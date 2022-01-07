#!/usr/bin/env python
__version__ = "2.1"
import argparse
import numpy as np
import vcf
import csv
import sys
import re
import math
import yaml
from collections import OrderedDict
from pathlib import Path
import subprocess
import os
import tempfile
import shutil
from collections import defaultdict
import configparser


def config_getter():
    this_dir = os.path.dirname(os.path.realpath(__file__))
    config_default = os.path.join(this_dir, "config.ini")
    config = configparser.ConfigParser(os.environ)
    config.read(config_default)
    hg38_ref = config.get("fasta", "hg38")
    hg19_ref = config.get("fasta", "hg19")
    hg38_anno = config.get("annotation", "hg38")
    hg19_anno = config.get("annotation", "hg19")
    yara_index = config.get("yara", "index")
    threads = config.get("yara", "threads")
    return {
        "hg38_ref": hg38_ref,
        "hg19_ref": hg19_ref,
        "hg38_anno": hg38_anno,
        "hg19_anno": hg19_anno,
        "yara_index": yara_index,
        "threads": threads,
    }


def is_insertion(ref, alt):
    return len(alt) > len(ref)


def is_deletion(ref, alt):
    return len(alt) < len(ref)


def simplify_indel_allele(ref, alt):
    while len(ref) > 0 and len(alt) > 0 and ref[-1] == alt[-1]:
        ref = ref[0:-1]
        alt = alt[0:-1]
    while len(ref) > 0 and len(alt) > 0 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
    return ref, alt


def parse_csq_format(vcf_reader):
    info_fields = vcf_reader.infos

    if info_fields["CSQ"] is None:
        sys.exit("Failed to extract format string from info description for tag (CSQ)")
    else:
        csq_header = info_fields["CSQ"]
        format_pattern = re.compile("Format: (.*)")
        match = format_pattern.search(csq_header.desc)
        return match.group(1)


def resolve_alleles(entry):
    alleles = {}
    if entry.is_indel:
        for alt in entry.ALT:
            alt = str(alt)
            if alt[0:1] != entry.REF[0:1]:
                alleles[alt] = alt
            elif alt[1:] == "":
                alleles[alt] = "-"
            else:
                alleles[alt] = alt[1:]
    elif "SVTYPE" in entry.INFO:
        if (
            entry.INFO["SVTYPE"] == "DEL" or entry.INFO["SVTYPE"] == "INS"
        ):  # PyVCF bug TYW
            for alt in entry.ALT:
                alt = str(alt)
                if alt[0:1] != entry.REF[0:1]:
                    alleles[alt] = alt
                elif alt[1:] == "":
                    alleles[alt] = "-"
                else:
                    alleles[alt] = alt[1:]

    else:
        for alt in entry.ALT:
            alt = str(alt)
            alleles[alt] = alt

    return alleles


def parse_csq_entries_for_allele(csq_entries, csq_format, csq_allele):
    csq_format_array = csq_format.split("|")

    transcripts = []
    for entry in csq_entries:
        values = entry.split("|")
        transcript = {}
        for key, value in zip(csq_format_array, values):
            transcript[key] = value
        if transcript["Allele"] == csq_allele:
            transcripts.append(transcript)

    return transcripts


def resolve_consequence(consequence_string):
    consequences = {
        consequence.lower() for consequence in consequence_string.split("&")
    }
    if "start_lost" in consequences:
        consequence = None
    elif "frameshift_variant" in consequences:
        consequence = "frameshift"
    elif "missense_variant" in consequences:
        consequence = "missense"
    elif "inframe_insertion" in consequences:
        consequence = "inframe_ins"
    elif "inframe_deletion" in consequences:
        consequence = "inframe_del"
    else:
        consequence = None
    return consequence


def calculate_coverage(ref, var):
    return ref + var


def calculate_vaf(ref, var):
    return (var / (calculate_coverage(ref, var) + 0.00001)) * 100


def convert_vcf(args_input=sys.argv[1:]):
    parser = argparse.ArgumentParser("convert_vcf")
    parser.add_argument(
        "input_file",
        type=argparse.FileType("r"),
        help="input VCF",
    )
    parser.add_argument(
        "output_file", type=argparse.FileType("w"), help="output list of variants"
    )
    args = parser.parse_args(args_input)

    vcf_reader = vcf.Reader(args.input_file)
    if len(vcf_reader.samples) > 1:
        sys.exit("ERROR: VCF file contains more than one sample")
    output_headers = [
        "chromosome_name",
        "start",
        "stop",
        "reference",
        "variant",
        "gene_name",
        "transcript_name",
        "amino_acid_change",
        "ensembl_gene_id",
        "wildtype_amino_acid_sequence",
        "downstream_amino_acid_sequence",
        "variant_type",
        "protein_position",
        "index",
    ]

    tsv_writer = csv.DictWriter(
        args.output_file, delimiter="\t", fieldnames=output_headers
    )
    tsv_writer.writeheader()
    csq_format = parse_csq_format(vcf_reader)
    transcript_count = {}
    for entry in vcf_reader:
        chromosome = entry.CHROM
        start = entry.affected_start
        stop = entry.affected_end
        reference = entry.REF
        alts = entry.ALT

        if len(vcf_reader.samples) == 1:
            genotype = entry.genotype(vcf_reader.samples[0])
            if genotype.gt_type is None or genotype.gt_type == 0:
                # The genotype is uncalled or hom_ref
                continue

        alleles_dict = resolve_alleles(entry)
        for alt in alts:
            alt = str(alt)
            if entry.is_indel:
                if is_deletion(reference, alt):
                    bam_readcount_position = start + 1
                    (simplified_reference, simplified_alt) = simplify_indel_allele(
                        reference, alt
                    )
                    ref_base = reference[1:2]
                    var_base = "-" + simplified_reference
                elif is_insertion(reference, alt):
                    bam_readcount_position = start
                    (simplified_reference, simplified_alt) = simplify_indel_allele(
                        reference, alt
                    )
                    ref_base = reference
                    var_base = "+" + simplified_alt
                variant_type = "indels"
            else:
                bam_readcount_position = entry.POS
                variant_type = "snvs"
                ref_base = reference
                var_base = alt

            csq_allele = alleles_dict[alt]
            transcripts = parse_csq_entries_for_allele(
                entry.INFO["CSQ"], csq_format, csq_allele
            )
            for transcript in transcripts:
                transcript_name = transcript["Feature"]
                if transcript_name in transcript_count:
                    transcript_count[transcript_name] += 1
                else:
                    transcript_count[transcript_name] = 1
                consequence = resolve_consequence(transcript["Consequence"])
                if consequence is None:
                    continue
                elif consequence == "frameshift":
                    if transcript["DownstreamProtein"] == "":
                        continue
                    else:
                        amino_acid_change_position = transcript["Protein_position"]
                else:
                    amino_acid_change_position = (
                        transcript["Protein_position"] + transcript["Amino_acids"]
                    )
                gene_name = transcript["SYMBOL"]
                index = f"{gene_name}_{transcript_name}_{transcript_count[transcript_name]}.{consequence}.{amino_acid_change_position}"
                ensembl_gene_id = transcript["Gene"]
                output_row = {
                    "chromosome_name": entry.CHROM,
                    "start": entry.affected_start,
                    "stop": entry.affected_end,
                    "reference": entry.REF,
                    "variant": alt,
                    "gene_name": gene_name,
                    "transcript_name": transcript_name,
                    "amino_acid_change": transcript["Amino_acids"],
                    "ensembl_gene_id": ensembl_gene_id,
                    "wildtype_amino_acid_sequence": transcript["WildtypeProtein"],
                    "downstream_amino_acid_sequence": transcript["DownstreamProtein"],
                    "variant_type": consequence,
                    "protein_position": transcript["Protein_position"],
                    "index": index,
                }
                if transcript["Amino_acids"]:
                    output_row["amino_acid_change"] = transcript["Amino_acids"]
                else:
                    output_row["amino_acid_change"] = "NA"

                tsv_writer.writerow(output_row)

    args.input_file.close()
    args.output_file.close()


###################################################################
# generate fasta functions

csv.field_size_limit(sys.maxsize)


def position_out_of_bounds(position, sequence):
    return position > len(sequence) - 1


# This subroutine is a bit funky but it was designed that way to mirror
# distance_from_end to increase code readability from the caller's perspective
def distance_from_start(position, string):
    return position


def distance_from_end(position, string):
    return len(string) - 1 - position


def determine_peptide_sequence_length(
    full_wildtype_sequence_length, peptide_sequence_length, line
):
    actual_peptide_sequence_length = peptide_sequence_length

    # If the wildtype sequence is shorter than the desired peptide sequence
    # length we use the wildtype sequence length instead so that the extraction
    # algorithm below works correctly
    if full_wildtype_sequence_length < actual_peptide_sequence_length:
        actual_peptide_sequence_length = full_wildtype_sequence_length
        print(
            "Wildtype sequence length is shorter than desired peptide sequence length at position ({}, {}, {}). Using wildtype sequence length ({}) instead.".format(
                line["chromosome_name"],
                line["start"],
                line["stop"],
                actual_peptide_sequence_length,
            )
        )

    return actual_peptide_sequence_length


def determine_flanking_sequence_length(
    full_wildtype_sequence_length, peptide_sequence_length, line
):
    actual_peptide_sequence_length = determine_peptide_sequence_length(
        full_wildtype_sequence_length, peptide_sequence_length, line
    )
    if actual_peptide_sequence_length % 2 == 0:
        return (actual_peptide_sequence_length - 2) / 2
    else:
        return (actual_peptide_sequence_length - 1) / 2


def get_wildtype_subsequence(
    position,
    full_wildtype_sequence,
    wildtype_amino_acid_length,
    peptide_sequence_length,
    line,
):
    one_flanking_sequence_length = int(
        determine_flanking_sequence_length(
            len(full_wildtype_sequence), peptide_sequence_length, line
        )
    )
    peptide_sequence_length = min(
        2 * one_flanking_sequence_length + wildtype_amino_acid_length,
        len(full_wildtype_sequence),
    )

    # We want to extract a subset from full_wildtype_sequence that is
    # peptide_sequence_length long so that the position ends
    # up in the middle of the extracted sequence.
    # If the position is too far toward the beginning or end of
    # full_wildtype_sequence there aren't enough amino acids on one side
    # to achieve this.
    if (
        distance_from_start(position, full_wildtype_sequence)
        < one_flanking_sequence_length
    ):
        wildtype_subsequence = full_wildtype_sequence[:peptide_sequence_length]
        mutation_position = position
    elif (
        distance_from_end(position, full_wildtype_sequence)
        < one_flanking_sequence_length
    ):
        start_position = len(full_wildtype_sequence) - peptide_sequence_length
        wildtype_subsequence = full_wildtype_sequence[start_position:]
        mutation_position = (
            peptide_sequence_length
            - distance_from_end(position, full_wildtype_sequence)
            - 1
        )
    elif (
        distance_from_start(position, full_wildtype_sequence)
        >= one_flanking_sequence_length
        and distance_from_end(position, full_wildtype_sequence)
        >= one_flanking_sequence_length
    ):
        start_position = position - one_flanking_sequence_length
        end_position = start_position + peptide_sequence_length
        wildtype_subsequence = full_wildtype_sequence[start_position:end_position]
        mutation_position = one_flanking_sequence_length
    else:
        sys.exit(
            "ERROR: Something went wrong during the retrieval of the wildtype sequence at position(%s, %s, %s)"
            % line["chromsome_name"],
            line["start"],
            line["stop"],
        )

    return mutation_position, wildtype_subsequence


def get_frameshift_subsequences(
    position, full_wildtype_sequence, peptide_sequence_length, line
):
    one_flanking_sequence_length = determine_flanking_sequence_length(
        len(full_wildtype_sequence), peptide_sequence_length, line
    )
    if position < one_flanking_sequence_length:
        start_position = 0
    else:
        start_position = int(position - one_flanking_sequence_length)
    wildtype_subsequence_stop_position = int(position + one_flanking_sequence_length)
    mutation_subsequence_stop_position = int(position)
    wildtype_subsequence = full_wildtype_sequence[
        start_position:wildtype_subsequence_stop_position
    ]
    mutation_start_subsequence = full_wildtype_sequence[
        start_position:mutation_subsequence_stop_position
    ]
    return start_position, wildtype_subsequence, mutation_start_subsequence


def combine_conflicting_variants(codon_changes):
    codon = list(codon_changes[0].split("/")[0].lower())
    modified_positions = []
    for codon_change in codon_changes:
        (old_codon, new_codon) = codon_change.split("/")
        change_positions = [
            i for i in range(len(old_codon)) if old_codon[i] != new_codon[i]
        ]
        for position in change_positions:
            if position in modified_positions:
                print("Warning: position has already been modified")
            codon[position] = new_codon[position].lower()
            modified_positions.append(position)
    return translate("".join(codon))


def generate_fasta(args_input=sys.argv[1:]):
    parser = argparse.ArgumentParser("generate_fasta")
    parser.add_argument(
        "input_file",
        type=argparse.FileType("r"),
        help="input list of variants",
    )
    parser.add_argument(
        "peptide_sequence_length", type=int, help="length of the peptide sequence"
    )
    parser.add_argument(
        "epitope_length", type=int, help="length of subpeptides(epitopes) to predict"
    )
    parser.add_argument(
        "output_file", type=argparse.FileType("w"), help="output FASTA file"
    )
    parser.add_argument(
        "output_key_file", type=argparse.FileType("w"), help="output FASTA key file"
    )
    parser.add_argument(
        "downstream_sequence_length",
        type=int,
        help="Cap to limit the downstream sequence length for frameshifts when creating the fasta file.",
    )
    args = parser.parse_args(args_input)

    peptide_sequence_length = args.peptide_sequence_length
    tsvin = csv.DictReader(args.input_file, delimiter="\t")
    fasta_sequences = OrderedDict()
    for line in tsvin:
        variant_type = line["variant_type"]
        full_wildtype_sequence = line["wildtype_amino_acid_sequence"]
        if variant_type == "frameshift":
            if (
                line["amino_acid_change"] is not None
                and line["amino_acid_change"].split("/")[0] == "-"
            ):
                position = int(line["protein_position"].split("-", 1)[0])
            else:
                position = int(line["protein_position"].split("-", 1)[0]) - 1
        elif variant_type == "missense" or variant_type == "inframe_ins":
            wildtype_amino_acid, mutant_amino_acid = line["amino_acid_change"].split(
                "/"
            )
            # deal with stop codon
            if "*" in wildtype_amino_acid:
                wildtype_amino_acid = wildtype_amino_acid.split("*")[0]
            elif "X" in wildtype_amino_acid:
                wildtype_amino_acid = wildtype_amino_acid.split("X")[0]
            if "*" in mutant_amino_acid:
                mutant_amino_acid = mutant_amino_acid.split("*")[0]
                stop_codon_added = True
            elif "X" in mutant_amino_acid:
                mutant_amino_acid = mutant_amino_acid.split("X")[0]
                stop_codon_added = True
            else:
                stop_codon_added = False

            if wildtype_amino_acid == "-":
                position = int(line["protein_position"].split("-", 1)[0])
                wildtype_amino_acid_length = 0
            else:
                if "-" in line["protein_position"]:
                    position = int(line["protein_position"].split("-", 1)[0]) - 1
                    wildtype_amino_acid_length = len(wildtype_amino_acid)
                else:
                    position = int(line["protein_position"]) - 1
                    wildtype_amino_acid_length = len(wildtype_amino_acid)
        elif variant_type == "inframe_del":
            variant_type = "inframe_del"
            wildtype_amino_acid, mutant_amino_acid = line["amino_acid_change"].split(
                "/"
            )
            # deal with stop codon
            if "*" in wildtype_amino_acid:
                wildtype_amino_acid = wildtype_amino_acid.split("*")[0]
            elif "X" in wildtype_amino_acid:
                wildtype_amino_acid = wildtype_amino_acid.split("X")[0]
            if "*" in mutant_amino_acid:
                mutant_amino_acid = mutant_amino_acid.split("*")[0]
                stop_codon_added = True
            elif "X" in mutant_amino_acid:
                mutant_amino_acid = mutant_amino_acid.split("X")[0]
                stop_codon_added = True
            else:
                stop_codon_added = False
            position = int(line["protein_position"].split("-", 1)[0]) - 1
            wildtype_amino_acid_length = len(wildtype_amino_acid)
            if mutant_amino_acid == "-":
                mutant_amino_acid = ""
        else:
            continue

        if position_out_of_bounds(position, full_wildtype_sequence):
            continue

        if variant_type == "frameshift":
            (
                mutation_start_position,
                wildtype_subsequence,
                mutant_subsequence,
            ) = get_frameshift_subsequences(
                position, full_wildtype_sequence, peptide_sequence_length, line
            )
            downstream_sequence = line["downstream_amino_acid_sequence"]

            if (
                args.downstream_sequence_length != 0
                and len(downstream_sequence) > args.downstream_sequence_length
            ):
                downstream_sequence = downstream_sequence[
                    0 : args.downstream_sequence_length
                ]
            mutant_subsequence += downstream_sequence
        else:
            mutation_start_position, wildtype_subsequence = get_wildtype_subsequence(
                position,
                full_wildtype_sequence,
                wildtype_amino_acid_length,
                peptide_sequence_length,
                line,
            )
            mutation_end_position = mutation_start_position + wildtype_amino_acid_length
            if (
                wildtype_amino_acid != "-"
                and wildtype_amino_acid
                != wildtype_subsequence[mutation_start_position:mutation_end_position]
            ):
                if line["amino_acid_change"].split("/")[0].count("*") > 1:
                    print(
                        "Warning: Amino acid change is not sane - contains multiple stops. Skipping entry {}".format(
                            line["index"]
                        )
                    )
                    continue
                else:
                    print(
                        "Warning: There was a mismatch between the actual wildtype amino acid and the expected amino acid. Did you use the same reference build version for VEP that you used for creating the VCF?\n%s"
                        % line
                    )
                    continue
            if stop_codon_added:
                mutant_subsequence = (
                    wildtype_subsequence[:mutation_start_position] + mutant_amino_acid
                )
            else:
                mutant_subsequence = (
                    wildtype_subsequence[:mutation_start_position]
                    + mutant_amino_acid
                    + wildtype_subsequence[mutation_end_position:]
                )

        if "*" in wildtype_subsequence or "*" in mutant_subsequence:
            continue

        if "X" in wildtype_subsequence or "X" in mutant_subsequence:
            continue

        if "U" in wildtype_subsequence or "U" in mutant_subsequence:
            print(
                "Warning. Sequence contains unsupported amino acid U. Skipping entry {}".format(
                    line["index"]
                )
            )
            continue

        if mutant_subsequence in wildtype_subsequence:
            # This is not a novel peptide
            continue

        if (
            len(wildtype_subsequence) < args.epitope_length
            or len(mutant_subsequence) < args.epitope_length
        ):
            continue

        variant_id = line["index"]
        for designation, subsequence in zip(
            ["WT", "MT"], [wildtype_subsequence, mutant_subsequence]
        ):
            # for designation, subsequence in zip(['MT'], [mutant_subsequence]):
            key = f"{designation}.{variant_id}"
            if subsequence in fasta_sequences:
                fasta_sequences[subsequence].append(key)
            else:
                fasta_sequences[subsequence] = [key]

    count = 1
    for (subsequence, keys) in list(fasta_sequences.items()):
        args.output_file.writelines(">%s\n" % count)
        args.output_file.writelines("%s\n" % subsequence)
        yaml.dump({count: keys}, args.output_key_file, default_flow_style=False)
        count += 1

    args.input_file.close()
    args.output_file.close()
    args.output_key_file.close()


#############parse output###########

csv.field_size_limit(sys.maxsize)


def parse_input_tsv_file(input_tsv_file):
    tsv_reader = csv.DictReader(input_tsv_file, delimiter="\t")
    tsv_entries = {}
    for line in tsv_reader:
        tsv_entries[line["index"]] = line
    return tsv_entries


def min_match_count(peptide_length):
    return math.ceil(peptide_length / 2)


def determine_consecutive_matches_from_left(mt_epitope_seq, wt_epitope_seq):
    consecutive_matches = 0
    for a, b in zip(mt_epitope_seq, wt_epitope_seq):
        if a == b:
            consecutive_matches += 1
        else:
            break
    return consecutive_matches


def determine_consecutive_matches_from_right(mt_epitope_seq, wt_epitope_seq):
    consecutive_matches = 0
    for a, b in zip(reversed(mt_epitope_seq), reversed(wt_epitope_seq)):
        if a == b:
            consecutive_matches += 1
        else:
            break
    return consecutive_matches


def determine_total_matches(mt_epitope_seq, wt_epitope_seq):
    matches = 0
    for a, b in zip(mt_epitope_seq, wt_epitope_seq):
        if a == b:
            matches += 1
    return matches


def find_mutation_position(wt_epitope_seq, mt_epitope_seq):
    for i, (wt_aa, mt_aa) in enumerate(zip(wt_epitope_seq, mt_epitope_seq)):
        if wt_aa != mt_aa:
            return i + 1
    return 0


def match_wildtype_and_mutant_entry_for_missense(
    result, mt_position, wt_results, previous_result
):
    # The WT epitope at the same position is the match
    match_position = mt_position
    mt_epitope_seq = result["mt_epitope_seq"]
    wt_result = wt_results[match_position]
    wt_epitope_seq = wt_result["wt_epitope_seq"]
    total_matches = determine_total_matches(mt_epitope_seq, wt_epitope_seq)
    if total_matches >= min_match_count(int(result["peptide_length"])):
        result["wt_epitope_seq"] = wt_epitope_seq
        result["wt_scores"] = wt_result["wt_scores"]
    else:
        result["wt_epitope_seq"] = "NA"
        result["wt_scores"] = dict.fromkeys(list(result["mt_scores"].keys()), "NA")

    if mt_epitope_seq == wt_epitope_seq:
        result["mutation_position"] = "NA"
    else:
        if previous_result:
            previous_mutation_position = previous_result["mutation_position"]
            if previous_mutation_position == "NA":
                result["mutation_position"] = find_mutation_position(
                    wt_epitope_seq, mt_epitope_seq
                )
            elif previous_mutation_position > 0:
                result["mutation_position"] = previous_mutation_position - 1
            else:
                result["mutation_position"] = 0
        else:
            result["mutation_position"] = find_mutation_position(
                wt_epitope_seq, mt_epitope_seq
            )


def match_wildtype_and_mutant_entry_for_frameshift(
    result, mt_position, wt_results, previous_result
):
    # The WT epitope at the same position is the match
    match_position = mt_position
    # Since the MT sequence is longer than the WT sequence, not all MT epitopes have a match
    if match_position not in wt_results:
        result["wt_epitope_seq"] = "NA"
        result["wt_scores"] = dict.fromkeys(list(result["mt_scores"].keys()), "NA")
        if previous_result["mutation_position"] == "NA":
            result["mutation_position"] = "NA"
        elif previous_result["mutation_position"] > 0:
            result["mutation_position"] = previous_result["mutation_position"] - 1
        else:
            result["mutation_position"] = 0
        return

    mt_epitope_seq = result["mt_epitope_seq"]
    wt_result = wt_results[match_position]
    wt_epitope_seq = wt_result["wt_epitope_seq"]
    if mt_epitope_seq == wt_epitope_seq:
        # The MT epitope does not overlap the frameshift mutation
        result["wt_epitope_seq"] = wt_result["wt_epitope_seq"]
        result["wt_scores"] = wt_result["wt_scores"]
        result["mutation_position"] = "NA"
    else:
        # Determine how many amino acids are the same between the MT epitope and its matching WT epitope
        total_matches = determine_total_matches(mt_epitope_seq, wt_epitope_seq)
        if total_matches >= min_match_count(int(result["peptide_length"])):
            # The minimum amino acid match count is met
            result["wt_epitope_seq"] = wt_result["wt_epitope_seq"]
            result["wt_scores"] = wt_result["wt_scores"]
        else:
            # The minimum amino acid match count is not met
            # Even though there is a matching WT epitope there are not enough overlapping amino acids
            # We don't include the matching WT epitope in the output
            result["wt_epitope_seq"] = "NA"
            result["wt_scores"] = dict.fromkeys(list(result["mt_scores"].keys()), "NA")
        mutation_position = find_mutation_position(wt_epitope_seq, mt_epitope_seq)
        # print('mu_pos', mutation_position)
        if mutation_position == 1 and int(previous_result["mutation_position"]) <= 1:
            # The true mutation position is to the left of the current MT eptiope
            mutation_position = 0
        result["mutation_position"] = mutation_position


def match_wildtype_and_mutant_entry_for_inframe_indel(
    result,
    mt_position,
    wt_results,
    previous_result,
    iedb_results_for_wt_iedb_result_key,
):
    # If the previous WT epitope was matched "from the right" we can just use that position to infer the mutation position and match direction
    if previous_result is not None and previous_result["match_direction"] == "right":
        best_match_position = previous_result["wt_epitope_position"] + 1
        result["wt_epitope_position"] = best_match_position
        result["match_direction"] = "right"
        if previous_result["mutation_position"] > 0:
            result["mutation_position"] = previous_result["mutation_position"] - 1
        else:
            result["mutation_position"] = 0

        # We need to ensure that the matched WT eptiope has enough overlapping amino acids with the MT epitope
        best_match_wt_result = wt_results[str(best_match_position)]
        total_matches = determine_total_matches(
            result["mt_epitope_seq"], best_match_wt_result["wt_epitope_seq"]
        )
        if total_matches and total_matches >= min_match_count(
            int(result["peptide_length"])
        ):
            # The minimum amino acid match count is met
            result["wt_epitope_seq"] = best_match_wt_result["wt_epitope_seq"]
            result["wt_scores"] = best_match_wt_result["wt_scores"]
        else:
            # The minimum amino acid match count is not met
            # Even though there is a matching WT epitope there are not enough overlapping amino acids
            # We don't include the matching WT epitope in the output
            result["wt_epitope_seq"] = "NA"
            result["wt_scores"] = dict.fromkeys(list(result["mt_scores"].keys()), "NA")

        return

    # In all other cases the WT epitope at the same position is used as the baseline match
    baseline_best_match_position = mt_position

    # For an inframe insertion the MT sequence is longer than the WT sequence
    # In this case not all MT epitopes might have a baseline match
    if baseline_best_match_position not in wt_results:
        result["wt_epitope_seq"] = "NA"
        result["wt_scores"] = dict.fromkeys(list(result["mt_scores"].keys()), "NA")
        # We then infer the mutation position and match direction from the previous MT epitope
        result["match_direction"] = previous_result["match_direction"]
        if previous_result["mutation_position"] > 0:
            result["mutation_position"] = previous_result["mutation_position"] - 1
        else:
            result["mutation_position"] = 0
        return

    mt_epitope_seq = result["mt_epitope_seq"]
    baseline_best_match_wt_result = wt_results[baseline_best_match_position]
    baseline_best_match_wt_epitope_seq = baseline_best_match_wt_result["wt_epitope_seq"]
    # The MT epitope does not overlap the indel mutation
    if baseline_best_match_wt_epitope_seq == mt_epitope_seq:
        result["wt_epitope_seq"] = baseline_best_match_wt_result["wt_epitope_seq"]
        result["wt_scores"] = baseline_best_match_wt_result["wt_scores"]
        result["wt_epitope_position"] = int(baseline_best_match_position)
        result["mutation_position"] = "NA"
        result["match_direction"] = "left"

    # If there is no previous result or the previous WT epitope was matched "from the left" we start by comparing to the baseline match
    if previous_result is None or previous_result["match_direction"] == "left":
        best_match_count = determine_consecutive_matches_from_left(
            mt_epitope_seq, baseline_best_match_wt_epitope_seq
        )
        # The alternate best match candidate "from the right" is inferred from the baseline best match position and the indel length
        if result["variant_type"] == "inframe_ins":
            insertion_length = len(
                list(iedb_results_for_wt_iedb_result_key.keys())
            ) - len(list(wt_results.keys()))
            alternate_best_match_position = (
                int(baseline_best_match_position) - insertion_length
            )
        elif result["variant_type"] == "inframe_del":
            deletion_length = len(list(wt_results.keys())) - len(
                list(iedb_results_for_wt_iedb_result_key.keys())
            )
            alternate_best_match_position = (
                int(baseline_best_match_position) + deletion_length
            )
        if alternate_best_match_position > 0:
            alternate_best_match_wt_result = wt_results[
                str(alternate_best_match_position)
            ]
            alternate_best_match_wt_epitope_seq = alternate_best_match_wt_result[
                "wt_epitope_seq"
            ]
            consecutive_matches_from_right = determine_consecutive_matches_from_right(
                mt_epitope_seq, alternate_best_match_wt_epitope_seq
            )
            # We then check if the alternate best match epitope has more matching amino acids than the baseline best match epitope
            # If it does, we pick it as the best match
            if consecutive_matches_from_right > best_match_count:
                match_direction = "right"
                best_match_position = alternate_best_match_position
                best_match_wt_result = alternate_best_match_wt_result
            else:
                match_direction = "left"
                best_match_position = baseline_best_match_position
                best_match_wt_result = baseline_best_match_wt_result
        else:
            match_direction = "left"
            best_match_position = baseline_best_match_position
            best_match_wt_result = baseline_best_match_wt_result

        # Now that we have found the matching WT epitope we still need to ensure that it has enough overlapping amino acids
        total_matches = determine_total_matches(
            mt_epitope_seq, best_match_wt_result["wt_epitope_seq"]
        )
        if total_matches and total_matches >= min_match_count(
            int(result["peptide_length"])
        ):
            # The minimum amino acid match count is met
            result["wt_epitope_seq"] = best_match_wt_result["wt_epitope_seq"]
            result["wt_scores"] = best_match_wt_result["wt_scores"]
        else:
            # The minimum amino acid match count is not met
            # Even though there is a matching WT epitope there are not enough overlapping amino acids
            # We don't include the matching WT epitope in the output
            result["wt_epitope_seq"] = "NA"
            result["wt_scores"] = dict.fromkeys(list(result["mt_scores"].keys()), "NA")

        result["mutation_position"] = find_mutation_position(
            baseline_best_match_wt_epitope_seq, mt_epitope_seq
        )
        result["match_direction"] = match_direction
        result["wt_epitope_position"] = best_match_position


def match_wildtype_and_mutant_entries(iedb_results, wt_iedb_results):
    new_iedb_results = {}

    for key in sorted(list(iedb_results.keys()), key=lambda x: int(x.split("|")[-1])):
        result = iedb_results[key]
        (wt_iedb_result_key, mt_position) = key.split("|", 1)

        previous_mt_position = str(int(mt_position) - 1)
        previous_key = "|".join([wt_iedb_result_key, previous_mt_position])

        if previous_key in iedb_results:
            previous_result = iedb_results[previous_key]
        else:
            previous_result = None
        wt_results = wt_iedb_results[wt_iedb_result_key]

        if result["variant_type"] == "missense":
            match_wildtype_and_mutant_entry_for_missense(
                result, mt_position, wt_results, previous_result
            )
        elif result["variant_type"] == "frameshift":
            match_wildtype_and_mutant_entry_for_frameshift(
                result, mt_position, wt_results, previous_result
            )
        elif (
            result["variant_type"] == "inframe_ins"
            or result["variant_type"] == "inframe_del"
        ):
            iedb_results_for_wt_iedb_result_key = {
                key: value
                for key, value in list(iedb_results.items())
                if key.startswith(wt_iedb_result_key)
            }
            match_wildtype_and_mutant_entry_for_inframe_indel(
                result,
                mt_position,
                wt_results,
                previous_result,
                iedb_results_for_wt_iedb_result_key,
            )

        mt_scores = result["mt_scores"]
        if min(mt_scores.values()) < 1000.0:
            new_iedb_results[key] = result
    return new_iedb_results


def parse_iedb_file(input_iedb_files, tsv_entries, key_file):
    protein_identifiers_from_label = yaml.safe_load(key_file)
    # print(protein_identifiers_from_label)
    iedb_results = {}
    wt_iedb_results = {}
    for input_iedb_file in input_iedb_files:
        iedb_tsv_reader = csv.DictReader(input_iedb_file, delimiter="\t")
        (method, remainder) = os.path.basename(input_iedb_file.name).split(".", 1)
        for line in iedb_tsv_reader:
            protein_label = int(line["seq_num"])
            if "core_peptide" in line and int(line["end"]) - int(line["start"]) == 8:
                # Start and end refer to the position of the core peptide
                # Infer the (start) position of the peptide from the positions of the core peptide
                position = str(
                    int(line["start"]) - line["peptide"].find(line["core_peptide"])
                )
            else:
                position = line["start"]
            epitope = line["peptide"]
            score = line["ic50"]
            allele = line["allele"]
            peptide_length = len(epitope)

            if protein_identifiers_from_label[protein_label] is not None:
                protein_identifiers = protein_identifiers_from_label[protein_label]

            for protein_identifier in protein_identifiers:
                (protein_type, tsv_index) = protein_identifier.split(".", 1)
                if protein_type == "MT":
                    if float(score) > 0.0:
                        tsv_entry = tsv_entries[tsv_index]
                        key = f"{tsv_index}|{position}"
                        if key not in iedb_results:
                            iedb_results[key] = {}
                            iedb_results[key]["mt_scores"] = {}
                            iedb_results[key]["mt_epitope_seq"] = epitope
                            iedb_results[key]["gene_name"] = tsv_entry["gene_name"]
                            iedb_results[key]["amino_acid_change"] = tsv_entry[
                                "amino_acid_change"
                            ]
                            iedb_results[key]["variant_type"] = tsv_entry[
                                "variant_type"
                            ]
                            iedb_results[key]["position"] = position
                            iedb_results[key]["tsv_index"] = tsv_index
                            iedb_results[key]["allele"] = allele
                            iedb_results[key]["peptide_length"] = peptide_length
                        iedb_results[key]["mt_scores"][method] = float(score)

                if protein_type == "WT":
                    if tsv_index not in wt_iedb_results:
                        wt_iedb_results[tsv_index] = {}
                    if position not in wt_iedb_results[tsv_index]:
                        wt_iedb_results[tsv_index][position] = {}
                        wt_iedb_results[tsv_index][position]["wt_scores"] = {}
                    wt_iedb_results[tsv_index][position]["wt_epitope_seq"] = epitope
                    wt_iedb_results[tsv_index][position]["wt_scores"][method] = float(
                        score
                    )

    return match_wildtype_and_mutant_entries(iedb_results, wt_iedb_results)


def add_summary_metrics(iedb_results):
    iedb_results_with_metrics = {}
    for key, value in list(iedb_results.items()):
        mt_scores = value["mt_scores"]
        best_mt_score = sys.maxsize
        for method, score in list(mt_scores.items()):
            if score < best_mt_score:
                best_mt_score = score
                best_mt_score_method = method
        value["best_mt_score"] = best_mt_score
        value["corresponding_wt_score"] = value["wt_scores"][best_mt_score_method]
        value["best_mt_score_method"] = best_mt_score_method
        value["median_mt_score"] = np.median(list(mt_scores.values()))
        wt_scores_with_value = [
            score for score in list(value["wt_scores"].values()) if score != "NA"
        ]
        if not wt_scores_with_value:
            value["median_wt_score"] = "NA"
        else:
            value["median_wt_score"] = np.median(wt_scores_with_value)
        iedb_results_with_metrics[key] = value

    return iedb_results_with_metrics


def pick_top_results(iedb_results, top_score_metric):
    score_at_position = {}
    for key, value in list(iedb_results.items()):
        (tsv_index, position) = key.split("|", 1)
        if tsv_index not in list(score_at_position.keys()):
            score_at_position[tsv_index] = {}
        if top_score_metric == "median":
            score_at_position[tsv_index][position] = value["median_mt_score"]
        elif top_score_metric == "lowest":
            score_at_position[tsv_index][position] = value["best_mt_score"]

    filtered_iedb_results = {}
    for tsv_index, value in list(score_at_position.items()):
        top_score = sys.maxsize
        for position, score in sorted(list(value.items()), key=lambda x: x[1]):
            top_score_key = f"{tsv_index}|{position}"
            if (
                iedb_results[top_score_key]["wt_epitope_seq"]
                != iedb_results[top_score_key]["mt_epitope_seq"]
            ):
                filtered_iedb_results[top_score_key] = iedb_results[top_score_key]
                break

    return filtered_iedb_results


def flatten_iedb_results(iedb_results):
    # transform the iedb_results dictionary into a two-dimensional list
    flattened_iedb_results = list(
        (
            value["gene_name"],
            value["amino_acid_change"],
            value["position"],
            value["mutation_position"],
            value["mt_scores"],
            value["wt_scores"],
            value["wt_epitope_seq"],
            value["mt_epitope_seq"],
            value["tsv_index"],
            value["allele"],
            value["peptide_length"],
            value["best_mt_score"],
            value["corresponding_wt_score"],
            value["best_mt_score_method"],
            value["median_mt_score"],
            value["median_wt_score"],
        )
        for value in list(iedb_results.values())
    )

    return flattened_iedb_results


def process_input_iedb_file(
    input_iedb_files, tsv_entries, key_file, top_result_per_mutation, top_score_metric
):
    iedb_results = parse_iedb_file(input_iedb_files, tsv_entries, key_file)
    iedb_results_with_metrics = add_summary_metrics(iedb_results)
    if top_result_per_mutation == True:
        filtered_iedb_results = pick_top_results(
            iedb_results_with_metrics, top_score_metric
        )
        flattened_iedb_results = flatten_iedb_results(filtered_iedb_results)
    else:
        flattened_iedb_results = flatten_iedb_results(iedb_results_with_metrics)

    return flattened_iedb_results


def output_headers(methods):
    # addtional headers
    method_dict = {
        "ann": "NetMHC",
        "netmhcpan": "NetMHCpan",
        "smm": "SMM",
        "pickpocket": "PickPocket",
    }
    headers = [
        "Chromosome",
        "Start",
        "Stop",
        "Reference",
        "Variant",
        "Transcript",
        "Ensembl Gene ID",
        "Variant Type",
        "Mutation",
        "Protein Position",
        "Gene Name",
        "HLA Allele",
        "Peptide Length",
        "Sub-peptide Position",
        "Mutation Position",
        "MT Epitope Seq",
        "WT Epitope Seq",
        "Best MT Score Method",
        "Best MT Score",
        "Corresponding WT Score",
        "Corresponding Fold Change",
        "Median MT Score",
        "Median WT Score",
        "Median Fold Change",
    ]

    for method in methods:
        pretty_method = method_dict[method]
        headers.append("%s WT Score" % pretty_method)
        headers.append("%s MT Score" % pretty_method)

    return headers


def determine_prediction_methods(input_iedb_files):
    methods = set()
    for input_iedb_file in input_iedb_files:
        (method, remainder) = os.path.basename(input_iedb_file.name).split(".", 1)
        methods.add(method)

    return sorted(list(methods))


def parse_output(args_input=sys.argv[1:]):
    parser = argparse.ArgumentParser("parse_output")
    parser.add_argument(
        "input_iedb_files",
        type=argparse.FileType("r"),
        nargs="+",
        help="Raw output file from IEDB",
    )
    parser.add_argument(
        "input_tsv_file", type=argparse.FileType("r"), help="Input list of variants"
    )
    parser.add_argument(
        "key_file", type=argparse.FileType("r"), help="Key file for lookup of FASTA IDs"
    )
    parser.add_argument("output_file", help="Parsed output file")
    parser.add_argument(
        "-t",
        "--top-result-per-mutation",
        action="store_true",
        default=False,
        help="Output top scoring candidate per allele-length per mutation. Default: False",
    )
    parser.add_argument(
        "-m",
        "--top-score-metric",
        choices=["lowest", "median"],
        default="lowest",
        help="The ic50 scoring metric to use when filtering for the top scoring results. "
        + "lowest: Best MT Score - lowest MT ic50 binding score of all chosen prediction methods. "
        + "median: Median MT Score All Methods - median MT ic50 binding score of all chosen prediction methods. "
        + "Default: lowest",
    )
    args = parser.parse_args(args_input)

    methods = determine_prediction_methods(args.input_iedb_files)
    tmp_output_file = args.output_file + ".tmp"
    tmp_output_filehandle = open(tmp_output_file, "w")
    tsv_writer = csv.DictWriter(
        tmp_output_filehandle, delimiter="\t", fieldnames=output_headers(methods)
    )
    tsv_writer.writeheader()

    tsv_entries = parse_input_tsv_file(args.input_tsv_file)

    # iedb_result   = parse_iedb_file(args.input_iedb_files, tsv_entries, args.key_file)
    iedb_results = process_input_iedb_file(
        args.input_iedb_files,
        tsv_entries,
        args.key_file,
        args.top_result_per_mutation,
        args.top_score_metric,
    )

    # TODO
    method_dict = {
        "ann": "NetMHC",
        "netmhcpan": "NetMHCpan",
        "smm": "SMM",
        "pickpocket": "PickPocket",
    }
    for (
        gene_name,
        variant_aa,
        position,
        mutation_position,
        mt_scores,
        wt_scores,
        wt_epitope_seq,
        mt_epitope_seq,
        tsv_index,
        allele,
        peptide_length,
        best_mt_score,
        corresponding_wt_score,
        best_mt_score_method,
        median_mt_score,
        median_wt_score,
    ) in iedb_results:
        tsv_entry = tsv_entries[tsv_index]
        if mt_epitope_seq != wt_epitope_seq:
            if wt_epitope_seq == "NA":
                corresponding_fold_change = "NA"
            else:
                corresponding_fold_change = "%.3f" % (
                    corresponding_wt_score / best_mt_score
                )
            if median_wt_score == "NA":
                median_fold_change = "NA"
            else:
                median_fold_change = "%.3f" % (median_wt_score / median_mt_score)
            row = {
                "Chromosome": tsv_entry["chromosome_name"],
                "Start": tsv_entry["start"],
                "Stop": tsv_entry["stop"],
                "Reference": tsv_entry["reference"],
                "Variant": tsv_entry["variant"],
                "Transcript": tsv_entry["transcript_name"],
                "Ensembl Gene ID": tsv_entry["ensembl_gene_id"],
                "Variant Type": tsv_entry["variant_type"],
                "Mutation": variant_aa,
                "Protein Position": tsv_entry["protein_position"],
                "Gene Name": gene_name,
                "HLA Allele": allele,
                "Peptide Length": peptide_length,
                "Sub-peptide Position": position,
                "Mutation Position": mutation_position,
                "MT Epitope Seq": mt_epitope_seq,
                "WT Epitope Seq": wt_epitope_seq,
                "Best MT Score Method": method_dict[best_mt_score_method],
                "Best MT Score": best_mt_score,
                "Corresponding WT Score": corresponding_wt_score,
                "Corresponding Fold Change": corresponding_fold_change,
                "Median MT Score": median_mt_score,
                "Median WT Score": median_wt_score,
                "Median Fold Change": median_fold_change,
            }
            for method in methods:
                pretty_method = method_dict[method]
                if method in wt_scores:
                    row["%s WT Score" % pretty_method] = wt_scores[method]
                else:
                    row["%s WT Score" % pretty_method] = "NA"
                if method in mt_scores:
                    row["%s MT Score" % pretty_method] = mt_scores[method]
                else:
                    row["%s MT Score" % pretty_method] = "NA"

            tsv_writer.writerow(row)

    tmp_output_filehandle.close()
    os.rename(tmp_output_file, args.output_file)
    for file_handle in args.input_iedb_files:
        file_handle.close()
    args.input_tsv_file.close()
    args.key_file.close()

    return args.output_file


#############################################################################
def run_cmd(cmd, msg=None):
    """"""
    status_message(cmd)
    if "," in msg:
        begin, finish = msg.split(",")
        status_message(begin)
    else:
        finish = msg
    try:
        subprocess.check_output(
            cmd,
            shell=True,
            stderr=subprocess.STDOUT,
        )
    except subprocess.CalledProcessError as err:
        error_msg = f"Error happend!: {err}, {err.output}"
    else:
        error_msg = ""
    if not error_msg:
        status_message(finish)
        return True
    else:
        status_message(error_msg)
        return False


def status_message(msg):
    print(msg)
    sys.stdout.flush()


def optitype_parser(result_tsv):
    """parse optitype results"""
    with open(result_tsv) as f:
        f.readline()
        l = f.readline().rstrip("\n").split("\t")
        hla_type_list = list(set(l[1:-2]))
        human_hla_list = [f"HLA-{j}" for j in hla_type_list]
        human_hla_list.sort()
        return ",".join(human_hla_list)


def sam2fastq(in_sam):
    """trim /1 and /2"""
    # out = tempfile.NamedTemporaryFile(delete=False)
    outfilename = os.path.splitext(os.path.basename(in_sam))[0]
    out = open(f"{outfilename}.fq", "w")
    with open(in_sam) as f:
        for line in f:
            if line.startswith("@"):
                continue
            l = line.rstrip("\n").split("\t")
            seq_id = l[0]
            seq = l[9]
            quality = l[10]
            out.write(f"@{seq_id}\n")
            out.write(f"{seq}\n")
            out.write("+\n")
            out.write(f"{quality}\n")
    out.close()
    return f"{outfilename}.fq"


def OptiType_wrapper(in_bam, config=config_getter()):
    """
    input: RNA-seq BAM file
    """
    try:
        path = subprocess.check_output(
            ["which", "OptiTypePipeline.py"], stderr=subprocess.STDOUT
        )
        path = path.decode("utf8").rstrip()
    except subprocess.CalledProcessError:
        sys.stderr.write(
            "Checking for {0} : ERROR - could not find {0}".format(
                "OptiTypePipeline.py"
            )
        )
        sys.stderr.write("Exiting.")
        sys.exit(0)

    hla_index = config["yara_index"]
    threads = config["threads"]

    optitype_folder = Path(path).resolve().parents[0]
    # print(optitype_folder)
    hla_ref = "{}/{}".format(optitype_folder, "data/hla_reference_rna.fasta")
    config = "{}/{}".format(optitype_folder, "config.ini")

    name = os.path.splitext(os.path.basename(in_bam))[0]
    hla_name = f"{name}.hla"

    bam2fastq = "picard -Xms1g -Xmx20g SamToFastq VALIDATION_STRINGENCY=SILENT I={0} F={1}_1.fq F2={1}_2.fq".format(
        in_bam, name
    )

    fastq_mapto_hla = (
        "yara_mapper -e 3 -t {0} {1} {2}_1.fq {2}_2.fq -o {2}.READS.bam".format(
            threads, hla_index, name
        )
    )

    hla_reads_filter = 'sambamba view -F "not unmapped and first_of_pair" {0}.READS.bam -f sam -o {1}_1.sam && sambamba view -F "not unmapped and second_of_pair" {0}.READS.bam -f sam -o {1}_2.sam'.format(
        name, hla_name
    )

    # results file: optitype/hla_{name}_result.tsv
    optitype_cmd = "OptiTypePipeline.py -i {0}_1.fq {0}_2.fq --rna -v -o optitype/ -c {1} -p hla_{2}".format(
        hla_name, config, name
    )

    flag = False
    ret1 = run_cmd(bam2fastq, "HLA typing step 1: FASTQ file generated!")
    if ret1:
        ret2 = run_cmd(
            fastq_mapto_hla, "HLA typing step 2: FASTQ map to HLA region finished!"
        )
        if ret2:
            ret3 = run_cmd(hla_reads_filter, "HLA typing step 3:  HLA SAM generated!")
            if ret3:
                hla_1_fq = sam2fastq(f"{hla_name}_1.sam")
                hla_2_fq = sam2fastq(f"{hla_name}_2.sam")
                if hla_1_fq and hla_2_fq:
                    status_message("HLA typing step 4:  HLA FASTQ generated!")
                    ret4 = run_cmd(
                        optitype_cmd,
                        "HLA typing step 5:  Optitype HLA typing finished!",
                    )
                    flag = True
                    os.remove(f"{name}_1.fq")
                    os.remove(f"{name}_2.fq")
                    os.remove(f"{name}.READS.bam")
                    os.remove(f"{hla_name}_1.fq")
                    os.remove(f"{hla_name}_2.fq")
                    os.remove(f"{hla_name}_1.sam")
                    os.remove(f"{hla_name}_2.sam")
                    return f"optitype/hla_{name}_result.tsv"


def allele_checker(in_allele, method):
    """
    input: in_alleles list, method = 'ann' or 'netmhcpan'
    output: valid alleles for specific MHC method
    """
    alleles_dict = {
        "ann": frozenset(
            [
                "HLA-A*01:01",
                "HLA-A*02:01",
                "HLA-A*02:02",
                "HLA-A*02:03",
                "HLA-A*02:06",
                "HLA-A*02:11",
                "HLA-A*02:12",
                "HLA-A*02:16",
                "HLA-A*02:17",
                "HLA-A*02:19",
                "HLA-A*02:50",
                "HLA-A*03:01",
                "HLA-A*11:01",
                "HLA-A*23:01",
                "HLA-A*24:02",
                "HLA-A*24:03",
                "HLA-A*25:01",
                "HLA-A*26:01",
                "HLA-A*26:02",
                "HLA-A*26:03",
                "HLA-A*29:02",
                "HLA-A*30:01",
                "HLA-A*30:02",
                "HLA-A*31:01",
                "HLA-A*32:01",
                "HLA-A*32:07",
                "HLA-A*32:15",
                "HLA-A*33:01",
                "HLA-A*66:01",
                "HLA-A*68:01",
                "HLA-A*68:02",
                "HLA-A*68:23",
                "HLA-A*69:01",
                "HLA-A*80:01",
                "HLA-B*07:02",
                "HLA-B*08:01",
                "HLA-B*08:02",
                "HLA-B*08:03",
                "HLA-B*14:02",
                "HLA-B*15:01",
                "HLA-B*15:02",
                "HLA-B*15:03",
                "HLA-B*15:09",
                "HLA-B*15:17",
                "HLA-B*18:01",
                "HLA-B*27:05",
                "HLA-B*27:20",
                "HLA-B*35:01",
                "HLA-B*35:03",
                "HLA-B*38:01",
                "HLA-B*39:01",
                "HLA-B*40:01",
                "HLA-B*40:02",
                "HLA-B*40:13",
                "HLA-B*42:01",
                "HLA-B*44:02",
                "HLA-B*44:03",
                "HLA-B*45:01",
                "HLA-B*46:01",
                "HLA-B*48:01",
                "HLA-B*51:01",
                "HLA-B*53:01",
                "HLA-B*54:01",
                "HLA-B*57:01",
                "HLA-B*58:01",
                "HLA-B*73:01",
                "HLA-B*83:01",
                "HLA-C*03:03",
                "HLA-C*04:01",
                "HLA-C*05:01",
                "HLA-C*06:02",
                "HLA-C*07:01",
                "HLA-C*07:02",
                "HLA-C*08:02",
                "HLA-C*12:03",
                "HLA-C*14:02",
                "HLA-C*15:02",
                "HLA-E*01:01",
            ]
        ),
        "netmhcpan": frozenset(
            [
                "HLA-A*01:01",
                "HLA-A*01:02",
                "HLA-A*01:03",
                "HLA-A*01:06",
                "HLA-A*01:07",
                "HLA-A*01:08",
                "HLA-A*01:09",
                "HLA-A*01:10",
                "HLA-A*01:12",
                "HLA-A*01:13",
                "HLA-A*01:14",
                "HLA-A*01:17",
                "HLA-A*01:19",
                "HLA-A*01:20",
                "HLA-A*01:21",
                "HLA-A*01:23",
                "HLA-A*01:24",
                "HLA-A*01:25",
                "HLA-A*01:26",
                "HLA-A*01:28",
                "HLA-A*01:29",
                "HLA-A*01:30",
                "HLA-A*01:32",
                "HLA-A*01:33",
                "HLA-A*01:35",
                "HLA-A*01:36",
                "HLA-A*01:37",
                "HLA-A*01:38",
                "HLA-A*01:39",
                "HLA-A*01:40",
                "HLA-A*01:41",
                "HLA-A*01:42",
                "HLA-A*01:43",
                "HLA-A*01:44",
                "HLA-A*01:45",
                "HLA-A*01:46",
                "HLA-A*01:47",
                "HLA-A*01:48",
                "HLA-A*01:49",
                "HLA-A*01:50",
                "HLA-A*01:51",
                "HLA-A*01:54",
                "HLA-A*01:55",
                "HLA-A*01:58",
                "HLA-A*01:59",
                "HLA-A*01:60",
                "HLA-A*01:61",
                "HLA-A*01:62",
                "HLA-A*01:63",
                "HLA-A*01:64",
                "HLA-A*01:65",
                "HLA-A*01:66",
                "HLA-A*02:01",
                "HLA-A*02:02",
                "HLA-A*02:03",
                "HLA-A*02:04",
                "HLA-A*02:05",
                "HLA-A*02:06",
                "HLA-A*02:07",
                "HLA-A*02:08",
                "HLA-A*02:09",
                "HLA-A*02:10",
                "HLA-A*02:101",
                "HLA-A*02:102",
                "HLA-A*02:103",
                "HLA-A*02:104",
                "HLA-A*02:105",
                "HLA-A*02:106",
                "HLA-A*02:107",
                "HLA-A*02:108",
                "HLA-A*02:109",
                "HLA-A*02:11",
                "HLA-A*02:110",
                "HLA-A*02:111",
                "HLA-A*02:112",
                "HLA-A*02:114",
                "HLA-A*02:115",
                "HLA-A*02:116",
                "HLA-A*02:117",
                "HLA-A*02:118",
                "HLA-A*02:119",
                "HLA-A*02:12",
                "HLA-A*02:120",
                "HLA-A*02:121",
                "HLA-A*02:122",
                "HLA-A*02:123",
                "HLA-A*02:124",
                "HLA-A*02:126",
                "HLA-A*02:127",
                "HLA-A*02:128",
                "HLA-A*02:129",
                "HLA-A*02:13",
                "HLA-A*02:130",
                "HLA-A*02:131",
                "HLA-A*02:132",
                "HLA-A*02:133",
                "HLA-A*02:134",
                "HLA-A*02:135",
                "HLA-A*02:136",
                "HLA-A*02:137",
                "HLA-A*02:138",
                "HLA-A*02:139",
                "HLA-A*02:14",
                "HLA-A*02:140",
                "HLA-A*02:141",
                "HLA-A*02:142",
                "HLA-A*02:143",
                "HLA-A*02:144",
                "HLA-A*02:145",
                "HLA-A*02:146",
                "HLA-A*02:147",
                "HLA-A*02:148",
                "HLA-A*02:149",
                "HLA-A*02:150",
                "HLA-A*02:151",
                "HLA-A*02:152",
                "HLA-A*02:153",
                "HLA-A*02:154",
                "HLA-A*02:155",
                "HLA-A*02:156",
                "HLA-A*02:157",
                "HLA-A*02:158",
                "HLA-A*02:159",
                "HLA-A*02:16",
                "HLA-A*02:160",
                "HLA-A*02:161",
                "HLA-A*02:162",
                "HLA-A*02:163",
                "HLA-A*02:164",
                "HLA-A*02:165",
                "HLA-A*02:166",
                "HLA-A*02:167",
                "HLA-A*02:168",
                "HLA-A*02:169",
                "HLA-A*02:17",
                "HLA-A*02:170",
                "HLA-A*02:171",
                "HLA-A*02:172",
                "HLA-A*02:173",
                "HLA-A*02:174",
                "HLA-A*02:175",
                "HLA-A*02:176",
                "HLA-A*02:177",
                "HLA-A*02:178",
                "HLA-A*02:179",
                "HLA-A*02:18",
                "HLA-A*02:180",
                "HLA-A*02:181",
                "HLA-A*02:182",
                "HLA-A*02:183",
                "HLA-A*02:184",
                "HLA-A*02:185",
                "HLA-A*02:186",
                "HLA-A*02:187",
                "HLA-A*02:188",
                "HLA-A*02:189",
                "HLA-A*02:19",
                "HLA-A*02:190",
                "HLA-A*02:191",
                "HLA-A*02:192",
                "HLA-A*02:193",
                "HLA-A*02:194",
                "HLA-A*02:195",
                "HLA-A*02:196",
                "HLA-A*02:197",
                "HLA-A*02:198",
                "HLA-A*02:199",
                "HLA-A*02:20",
                "HLA-A*02:200",
                "HLA-A*02:201",
                "HLA-A*02:202",
                "HLA-A*02:203",
                "HLA-A*02:204",
                "HLA-A*02:205",
                "HLA-A*02:206",
                "HLA-A*02:207",
                "HLA-A*02:208",
                "HLA-A*02:209",
                "HLA-A*02:21",
                "HLA-A*02:210",
                "HLA-A*02:211",
                "HLA-A*02:212",
                "HLA-A*02:213",
                "HLA-A*02:214",
                "HLA-A*02:215",
                "HLA-A*02:216",
                "HLA-A*02:217",
                "HLA-A*02:218",
                "HLA-A*02:219",
                "HLA-A*02:22",
                "HLA-A*02:220",
                "HLA-A*02:221",
                "HLA-A*02:224",
                "HLA-A*02:228",
                "HLA-A*02:229",
                "HLA-A*02:230",
                "HLA-A*02:231",
                "HLA-A*02:232",
                "HLA-A*02:233",
                "HLA-A*02:234",
                "HLA-A*02:235",
                "HLA-A*02:236",
                "HLA-A*02:237",
                "HLA-A*02:238",
                "HLA-A*02:239",
                "HLA-A*02:24",
                "HLA-A*02:240",
                "HLA-A*02:241",
                "HLA-A*02:242",
                "HLA-A*02:243",
                "HLA-A*02:244",
                "HLA-A*02:245",
                "HLA-A*02:246",
                "HLA-A*02:247",
                "HLA-A*02:248",
                "HLA-A*02:249",
                "HLA-A*02:25",
                "HLA-A*02:251",
                "HLA-A*02:252",
                "HLA-A*02:253",
                "HLA-A*02:254",
                "HLA-A*02:255",
                "HLA-A*02:256",
                "HLA-A*02:257",
                "HLA-A*02:258",
                "HLA-A*02:259",
                "HLA-A*02:26",
                "HLA-A*02:260",
                "HLA-A*02:261",
                "HLA-A*02:262",
                "HLA-A*02:263",
                "HLA-A*02:264",
                "HLA-A*02:265",
                "HLA-A*02:266",
                "HLA-A*02:27",
                "HLA-A*02:28",
                "HLA-A*02:29",
                "HLA-A*02:30",
                "HLA-A*02:31",
                "HLA-A*02:33",
                "HLA-A*02:34",
                "HLA-A*02:35",
                "HLA-A*02:36",
                "HLA-A*02:37",
                "HLA-A*02:38",
                "HLA-A*02:39",
                "HLA-A*02:40",
                "HLA-A*02:41",
                "HLA-A*02:42",
                "HLA-A*02:44",
                "HLA-A*02:45",
                "HLA-A*02:46",
                "HLA-A*02:47",
                "HLA-A*02:48",
                "HLA-A*02:49",
                "HLA-A*02:50",
                "HLA-A*02:51",
                "HLA-A*02:52",
                "HLA-A*02:54",
                "HLA-A*02:55",
                "HLA-A*02:56",
                "HLA-A*02:57",
                "HLA-A*02:58",
                "HLA-A*02:59",
                "HLA-A*02:60",
                "HLA-A*02:61",
                "HLA-A*02:62",
                "HLA-A*02:63",
                "HLA-A*02:64",
                "HLA-A*02:65",
                "HLA-A*02:66",
                "HLA-A*02:67",
                "HLA-A*02:68",
                "HLA-A*02:69",
                "HLA-A*02:70",
                "HLA-A*02:71",
                "HLA-A*02:72",
                "HLA-A*02:73",
                "HLA-A*02:74",
                "HLA-A*02:75",
                "HLA-A*02:76",
                "HLA-A*02:77",
                "HLA-A*02:78",
                "HLA-A*02:79",
                "HLA-A*02:80",
                "HLA-A*02:81",
                "HLA-A*02:84",
                "HLA-A*02:85",
                "HLA-A*02:86",
                "HLA-A*02:87",
                "HLA-A*02:89",
                "HLA-A*02:90",
                "HLA-A*02:91",
                "HLA-A*02:92",
                "HLA-A*02:93",
                "HLA-A*02:95",
                "HLA-A*02:96",
                "HLA-A*02:97",
                "HLA-A*02:99",
                "HLA-A*03:01",
                "HLA-A*03:02",
                "HLA-A*03:04",
                "HLA-A*03:05",
                "HLA-A*03:06",
                "HLA-A*03:07",
                "HLA-A*03:08",
                "HLA-A*03:09",
                "HLA-A*03:10",
                "HLA-A*03:12",
                "HLA-A*03:13",
                "HLA-A*03:14",
                "HLA-A*03:15",
                "HLA-A*03:16",
                "HLA-A*03:17",
                "HLA-A*03:18",
                "HLA-A*03:19",
                "HLA-A*03:20",
                "HLA-A*03:22",
                "HLA-A*03:23",
                "HLA-A*03:24",
                "HLA-A*03:25",
                "HLA-A*03:26",
                "HLA-A*03:27",
                "HLA-A*03:28",
                "HLA-A*03:29",
                "HLA-A*03:30",
                "HLA-A*03:31",
                "HLA-A*03:32",
                "HLA-A*03:33",
                "HLA-A*03:34",
                "HLA-A*03:35",
                "HLA-A*03:37",
                "HLA-A*03:38",
                "HLA-A*03:39",
                "HLA-A*03:40",
                "HLA-A*03:41",
                "HLA-A*03:42",
                "HLA-A*03:43",
                "HLA-A*03:44",
                "HLA-A*03:45",
                "HLA-A*03:46",
                "HLA-A*03:47",
                "HLA-A*03:48",
                "HLA-A*03:49",
                "HLA-A*03:50",
                "HLA-A*03:51",
                "HLA-A*03:52",
                "HLA-A*03:53",
                "HLA-A*03:54",
                "HLA-A*03:55",
                "HLA-A*03:56",
                "HLA-A*03:57",
                "HLA-A*03:58",
                "HLA-A*03:59",
                "HLA-A*03:60",
                "HLA-A*03:61",
                "HLA-A*03:62",
                "HLA-A*03:63",
                "HLA-A*03:64",
                "HLA-A*03:65",
                "HLA-A*03:66",
                "HLA-A*03:67",
                "HLA-A*03:70",
                "HLA-A*03:71",
                "HLA-A*03:72",
                "HLA-A*03:73",
                "HLA-A*03:74",
                "HLA-A*03:75",
                "HLA-A*03:76",
                "HLA-A*03:77",
                "HLA-A*03:78",
                "HLA-A*03:79",
                "HLA-A*03:80",
                "HLA-A*03:81",
                "HLA-A*03:82",
                "HLA-A*11:01",
                "HLA-A*11:02",
                "HLA-A*11:03",
                "HLA-A*11:04",
                "HLA-A*11:05",
                "HLA-A*11:06",
                "HLA-A*11:07",
                "HLA-A*11:08",
                "HLA-A*11:09",
                "HLA-A*11:10",
                "HLA-A*11:11",
                "HLA-A*11:12",
                "HLA-A*11:13",
                "HLA-A*11:14",
                "HLA-A*11:15",
                "HLA-A*11:16",
                "HLA-A*11:17",
                "HLA-A*11:18",
                "HLA-A*11:19",
                "HLA-A*11:20",
                "HLA-A*11:22",
                "HLA-A*11:23",
                "HLA-A*11:24",
                "HLA-A*11:25",
                "HLA-A*11:26",
                "HLA-A*11:27",
                "HLA-A*11:29",
                "HLA-A*11:30",
                "HLA-A*11:31",
                "HLA-A*11:32",
                "HLA-A*11:33",
                "HLA-A*11:34",
                "HLA-A*11:35",
                "HLA-A*11:36",
                "HLA-A*11:37",
                "HLA-A*11:38",
                "HLA-A*11:39",
                "HLA-A*11:40",
                "HLA-A*11:41",
                "HLA-A*11:42",
                "HLA-A*11:43",
                "HLA-A*11:44",
                "HLA-A*11:45",
                "HLA-A*11:46",
                "HLA-A*11:47",
                "HLA-A*11:48",
                "HLA-A*11:49",
                "HLA-A*11:51",
                "HLA-A*11:53",
                "HLA-A*11:54",
                "HLA-A*11:55",
                "HLA-A*11:56",
                "HLA-A*11:57",
                "HLA-A*11:58",
                "HLA-A*11:59",
                "HLA-A*11:60",
                "HLA-A*11:61",
                "HLA-A*11:62",
                "HLA-A*11:63",
                "HLA-A*11:64",
                "HLA-A*23:01",
                "HLA-A*23:02",
                "HLA-A*23:03",
                "HLA-A*23:04",
                "HLA-A*23:05",
                "HLA-A*23:06",
                "HLA-A*23:09",
                "HLA-A*23:10",
                "HLA-A*23:12",
                "HLA-A*23:13",
                "HLA-A*23:14",
                "HLA-A*23:15",
                "HLA-A*23:16",
                "HLA-A*23:17",
                "HLA-A*23:18",
                "HLA-A*23:20",
                "HLA-A*23:21",
                "HLA-A*23:22",
                "HLA-A*23:23",
                "HLA-A*23:24",
                "HLA-A*23:25",
                "HLA-A*23:26",
                "HLA-A*24:02",
                "HLA-A*24:03",
                "HLA-A*24:04",
                "HLA-A*24:05",
                "HLA-A*24:06",
                "HLA-A*24:07",
                "HLA-A*24:08",
                "HLA-A*24:10",
                "HLA-A*24:100",
                "HLA-A*24:101",
                "HLA-A*24:102",
                "HLA-A*24:103",
                "HLA-A*24:104",
                "HLA-A*24:105",
                "HLA-A*24:106",
                "HLA-A*24:107",
                "HLA-A*24:108",
                "HLA-A*24:109",
                "HLA-A*24:110",
                "HLA-A*24:111",
                "HLA-A*24:112",
                "HLA-A*24:113",
                "HLA-A*24:114",
                "HLA-A*24:115",
                "HLA-A*24:116",
                "HLA-A*24:117",
                "HLA-A*24:118",
                "HLA-A*24:119",
                "HLA-A*24:120",
                "HLA-A*24:121",
                "HLA-A*24:122",
                "HLA-A*24:123",
                "HLA-A*24:124",
                "HLA-A*24:125",
                "HLA-A*24:126",
                "HLA-A*24:127",
                "HLA-A*24:128",
                "HLA-A*24:129",
                "HLA-A*24:13",
                "HLA-A*24:130",
                "HLA-A*24:131",
                "HLA-A*24:133",
                "HLA-A*24:134",
                "HLA-A*24:135",
                "HLA-A*24:136",
                "HLA-A*24:137",
                "HLA-A*24:138",
                "HLA-A*24:139",
                "HLA-A*24:14",
                "HLA-A*24:140",
                "HLA-A*24:141",
                "HLA-A*24:142",
                "HLA-A*24:143",
                "HLA-A*24:144",
                "HLA-A*24:15",
                "HLA-A*24:17",
                "HLA-A*24:18",
                "HLA-A*24:19",
                "HLA-A*24:20",
                "HLA-A*24:21",
                "HLA-A*24:22",
                "HLA-A*24:23",
                "HLA-A*24:24",
                "HLA-A*24:25",
                "HLA-A*24:26",
                "HLA-A*24:27",
                "HLA-A*24:28",
                "HLA-A*24:29",
                "HLA-A*24:30",
                "HLA-A*24:31",
                "HLA-A*24:32",
                "HLA-A*24:33",
                "HLA-A*24:34",
                "HLA-A*24:35",
                "HLA-A*24:37",
                "HLA-A*24:38",
                "HLA-A*24:39",
                "HLA-A*24:41",
                "HLA-A*24:42",
                "HLA-A*24:43",
                "HLA-A*24:44",
                "HLA-A*24:46",
                "HLA-A*24:47",
                "HLA-A*24:49",
                "HLA-A*24:50",
                "HLA-A*24:51",
                "HLA-A*24:52",
                "HLA-A*24:53",
                "HLA-A*24:54",
                "HLA-A*24:55",
                "HLA-A*24:56",
                "HLA-A*24:57",
                "HLA-A*24:58",
                "HLA-A*24:59",
                "HLA-A*24:61",
                "HLA-A*24:62",
                "HLA-A*24:63",
                "HLA-A*24:64",
                "HLA-A*24:66",
                "HLA-A*24:67",
                "HLA-A*24:68",
                "HLA-A*24:69",
                "HLA-A*24:70",
                "HLA-A*24:71",
                "HLA-A*24:72",
                "HLA-A*24:73",
                "HLA-A*24:74",
                "HLA-A*24:75",
                "HLA-A*24:76",
                "HLA-A*24:77",
                "HLA-A*24:78",
                "HLA-A*24:79",
                "HLA-A*24:80",
                "HLA-A*24:81",
                "HLA-A*24:82",
                "HLA-A*24:85",
                "HLA-A*24:87",
                "HLA-A*24:88",
                "HLA-A*24:89",
                "HLA-A*24:91",
                "HLA-A*24:92",
                "HLA-A*24:93",
                "HLA-A*24:94",
                "HLA-A*24:95",
                "HLA-A*24:96",
                "HLA-A*24:97",
                "HLA-A*24:98",
                "HLA-A*24:99",
                "HLA-A*25:01",
                "HLA-A*25:02",
                "HLA-A*25:03",
                "HLA-A*25:04",
                "HLA-A*25:05",
                "HLA-A*25:06",
                "HLA-A*25:07",
                "HLA-A*25:08",
                "HLA-A*25:09",
                "HLA-A*25:10",
                "HLA-A*25:11",
                "HLA-A*25:13",
                "HLA-A*26:01",
                "HLA-A*26:02",
                "HLA-A*26:03",
                "HLA-A*26:04",
                "HLA-A*26:05",
                "HLA-A*26:06",
                "HLA-A*26:07",
                "HLA-A*26:08",
                "HLA-A*26:09",
                "HLA-A*26:10",
                "HLA-A*26:12",
                "HLA-A*26:13",
                "HLA-A*26:14",
                "HLA-A*26:15",
                "HLA-A*26:16",
                "HLA-A*26:17",
                "HLA-A*26:18",
                "HLA-A*26:19",
                "HLA-A*26:20",
                "HLA-A*26:21",
                "HLA-A*26:22",
                "HLA-A*26:23",
                "HLA-A*26:24",
                "HLA-A*26:26",
                "HLA-A*26:27",
                "HLA-A*26:28",
                "HLA-A*26:29",
                "HLA-A*26:30",
                "HLA-A*26:31",
                "HLA-A*26:32",
                "HLA-A*26:33",
                "HLA-A*26:34",
                "HLA-A*26:35",
                "HLA-A*26:36",
                "HLA-A*26:37",
                "HLA-A*26:38",
                "HLA-A*26:39",
                "HLA-A*26:40",
                "HLA-A*26:41",
                "HLA-A*26:42",
                "HLA-A*26:43",
                "HLA-A*26:45",
                "HLA-A*26:46",
                "HLA-A*26:47",
                "HLA-A*26:48",
                "HLA-A*26:49",
                "HLA-A*26:50",
                "HLA-A*29:01",
                "HLA-A*29:02",
                "HLA-A*29:03",
                "HLA-A*29:04",
                "HLA-A*29:05",
                "HLA-A*29:06",
                "HLA-A*29:07",
                "HLA-A*29:09",
                "HLA-A*29:10",
                "HLA-A*29:11",
                "HLA-A*29:12",
                "HLA-A*29:13",
                "HLA-A*29:14",
                "HLA-A*29:15",
                "HLA-A*29:16",
                "HLA-A*29:17",
                "HLA-A*29:18",
                "HLA-A*29:19",
                "HLA-A*29:20",
                "HLA-A*29:21",
                "HLA-A*29:22",
                "HLA-A*30:01",
                "HLA-A*30:02",
                "HLA-A*30:03",
                "HLA-A*30:04",
                "HLA-A*30:06",
                "HLA-A*30:07",
                "HLA-A*30:08",
                "HLA-A*30:09",
                "HLA-A*30:10",
                "HLA-A*30:11",
                "HLA-A*30:12",
                "HLA-A*30:13",
                "HLA-A*30:15",
                "HLA-A*30:16",
                "HLA-A*30:17",
                "HLA-A*30:18",
                "HLA-A*30:19",
                "HLA-A*30:20",
                "HLA-A*30:22",
                "HLA-A*30:23",
                "HLA-A*30:24",
                "HLA-A*30:25",
                "HLA-A*30:26",
                "HLA-A*30:28",
                "HLA-A*30:29",
                "HLA-A*30:30",
                "HLA-A*30:31",
                "HLA-A*30:32",
                "HLA-A*30:33",
                "HLA-A*30:34",
                "HLA-A*30:35",
                "HLA-A*30:36",
                "HLA-A*30:37",
                "HLA-A*30:38",
                "HLA-A*30:39",
                "HLA-A*30:40",
                "HLA-A*30:41",
                "HLA-A*31:01",
                "HLA-A*31:02",
                "HLA-A*31:03",
                "HLA-A*31:04",
                "HLA-A*31:05",
                "HLA-A*31:06",
                "HLA-A*31:07",
                "HLA-A*31:08",
                "HLA-A*31:09",
                "HLA-A*31:10",
                "HLA-A*31:11",
                "HLA-A*31:12",
                "HLA-A*31:13",
                "HLA-A*31:15",
                "HLA-A*31:16",
                "HLA-A*31:17",
                "HLA-A*31:18",
                "HLA-A*31:19",
                "HLA-A*31:20",
                "HLA-A*31:21",
                "HLA-A*31:22",
                "HLA-A*31:23",
                "HLA-A*31:24",
                "HLA-A*31:25",
                "HLA-A*31:26",
                "HLA-A*31:27",
                "HLA-A*31:28",
                "HLA-A*31:29",
                "HLA-A*31:30",
                "HLA-A*31:31",
                "HLA-A*31:32",
                "HLA-A*31:33",
                "HLA-A*31:34",
                "HLA-A*31:35",
                "HLA-A*31:36",
                "HLA-A*31:37",
                "HLA-A*32:01",
                "HLA-A*32:02",
                "HLA-A*32:03",
                "HLA-A*32:04",
                "HLA-A*32:05",
                "HLA-A*32:06",
                "HLA-A*32:07",
                "HLA-A*32:08",
                "HLA-A*32:09",
                "HLA-A*32:10",
                "HLA-A*32:12",
                "HLA-A*32:13",
                "HLA-A*32:14",
                "HLA-A*32:15",
                "HLA-A*32:16",
                "HLA-A*32:17",
                "HLA-A*32:18",
                "HLA-A*32:20",
                "HLA-A*32:21",
                "HLA-A*32:22",
                "HLA-A*32:23",
                "HLA-A*32:24",
                "HLA-A*32:25",
                "HLA-A*33:01",
                "HLA-A*33:03",
                "HLA-A*33:04",
                "HLA-A*33:05",
                "HLA-A*33:06",
                "HLA-A*33:07",
                "HLA-A*33:08",
                "HLA-A*33:09",
                "HLA-A*33:10",
                "HLA-A*33:11",
                "HLA-A*33:12",
                "HLA-A*33:13",
                "HLA-A*33:14",
                "HLA-A*33:15",
                "HLA-A*33:16",
                "HLA-A*33:17",
                "HLA-A*33:18",
                "HLA-A*33:19",
                "HLA-A*33:20",
                "HLA-A*33:21",
                "HLA-A*33:22",
                "HLA-A*33:23",
                "HLA-A*33:24",
                "HLA-A*33:25",
                "HLA-A*33:26",
                "HLA-A*33:27",
                "HLA-A*33:28",
                "HLA-A*33:29",
                "HLA-A*33:30",
                "HLA-A*33:31",
                "HLA-A*34:01",
                "HLA-A*34:02",
                "HLA-A*34:03",
                "HLA-A*34:04",
                "HLA-A*34:05",
                "HLA-A*34:06",
                "HLA-A*34:07",
                "HLA-A*34:08",
                "HLA-A*36:01",
                "HLA-A*36:02",
                "HLA-A*36:03",
                "HLA-A*36:04",
                "HLA-A*36:05",
                "HLA-A*43:01",
                "HLA-A*66:01",
                "HLA-A*66:02",
                "HLA-A*66:03",
                "HLA-A*66:04",
                "HLA-A*66:05",
                "HLA-A*66:06",
                "HLA-A*66:07",
                "HLA-A*66:08",
                "HLA-A*66:09",
                "HLA-A*66:10",
                "HLA-A*66:11",
                "HLA-A*66:12",
                "HLA-A*66:13",
                "HLA-A*66:14",
                "HLA-A*66:15",
                "HLA-A*68:01",
                "HLA-A*68:02",
                "HLA-A*68:03",
                "HLA-A*68:04",
                "HLA-A*68:05",
                "HLA-A*68:06",
                "HLA-A*68:07",
                "HLA-A*68:08",
                "HLA-A*68:09",
                "HLA-A*68:10",
                "HLA-A*68:12",
                "HLA-A*68:13",
                "HLA-A*68:14",
                "HLA-A*68:15",
                "HLA-A*68:16",
                "HLA-A*68:17",
                "HLA-A*68:19",
                "HLA-A*68:20",
                "HLA-A*68:21",
                "HLA-A*68:22",
                "HLA-A*68:23",
                "HLA-A*68:24",
                "HLA-A*68:25",
                "HLA-A*68:26",
                "HLA-A*68:27",
                "HLA-A*68:28",
                "HLA-A*68:29",
                "HLA-A*68:30",
                "HLA-A*68:31",
                "HLA-A*68:32",
                "HLA-A*68:33",
                "HLA-A*68:34",
                "HLA-A*68:35",
                "HLA-A*68:36",
                "HLA-A*68:37",
                "HLA-A*68:38",
                "HLA-A*68:39",
                "HLA-A*68:40",
                "HLA-A*68:41",
                "HLA-A*68:42",
                "HLA-A*68:43",
                "HLA-A*68:44",
                "HLA-A*68:45",
                "HLA-A*68:46",
                "HLA-A*68:47",
                "HLA-A*68:48",
                "HLA-A*68:50",
                "HLA-A*68:51",
                "HLA-A*68:52",
                "HLA-A*68:53",
                "HLA-A*68:54",
                "HLA-A*69:01",
                "HLA-A*74:01",
                "HLA-A*74:02",
                "HLA-A*74:03",
                "HLA-A*74:04",
                "HLA-A*74:05",
                "HLA-A*74:06",
                "HLA-A*74:07",
                "HLA-A*74:08",
                "HLA-A*74:09",
                "HLA-A*74:10",
                "HLA-A*74:11",
                "HLA-A*74:13",
                "HLA-A*80:01",
                "HLA-A*80:02",
                "HLA-B*07:02",
                "HLA-B*07:03",
                "HLA-B*07:04",
                "HLA-B*07:05",
                "HLA-B*07:06",
                "HLA-B*07:07",
                "HLA-B*07:08",
                "HLA-B*07:09",
                "HLA-B*07:10",
                "HLA-B*07:100",
                "HLA-B*07:101",
                "HLA-B*07:102",
                "HLA-B*07:103",
                "HLA-B*07:104",
                "HLA-B*07:105",
                "HLA-B*07:106",
                "HLA-B*07:107",
                "HLA-B*07:108",
                "HLA-B*07:109",
                "HLA-B*07:11",
                "HLA-B*07:110",
                "HLA-B*07:112",
                "HLA-B*07:113",
                "HLA-B*07:114",
                "HLA-B*07:115",
                "HLA-B*07:12",
                "HLA-B*07:13",
                "HLA-B*07:14",
                "HLA-B*07:15",
                "HLA-B*07:16",
                "HLA-B*07:17",
                "HLA-B*07:18",
                "HLA-B*07:19",
                "HLA-B*07:20",
                "HLA-B*07:21",
                "HLA-B*07:22",
                "HLA-B*07:23",
                "HLA-B*07:24",
                "HLA-B*07:25",
                "HLA-B*07:26",
                "HLA-B*07:27",
                "HLA-B*07:28",
                "HLA-B*07:29",
                "HLA-B*07:30",
                "HLA-B*07:31",
                "HLA-B*07:32",
                "HLA-B*07:33",
                "HLA-B*07:34",
                "HLA-B*07:35",
                "HLA-B*07:36",
                "HLA-B*07:37",
                "HLA-B*07:38",
                "HLA-B*07:39",
                "HLA-B*07:40",
                "HLA-B*07:41",
                "HLA-B*07:42",
                "HLA-B*07:43",
                "HLA-B*07:44",
                "HLA-B*07:45",
                "HLA-B*07:46",
                "HLA-B*07:47",
                "HLA-B*07:48",
                "HLA-B*07:50",
                "HLA-B*07:51",
                "HLA-B*07:52",
                "HLA-B*07:53",
                "HLA-B*07:54",
                "HLA-B*07:55",
                "HLA-B*07:56",
                "HLA-B*07:57",
                "HLA-B*07:58",
                "HLA-B*07:59",
                "HLA-B*07:60",
                "HLA-B*07:61",
                "HLA-B*07:62",
                "HLA-B*07:63",
                "HLA-B*07:64",
                "HLA-B*07:65",
                "HLA-B*07:66",
                "HLA-B*07:68",
                "HLA-B*07:69",
                "HLA-B*07:70",
                "HLA-B*07:71",
                "HLA-B*07:72",
                "HLA-B*07:73",
                "HLA-B*07:74",
                "HLA-B*07:75",
                "HLA-B*07:76",
                "HLA-B*07:77",
                "HLA-B*07:78",
                "HLA-B*07:79",
                "HLA-B*07:80",
                "HLA-B*07:81",
                "HLA-B*07:82",
                "HLA-B*07:83",
                "HLA-B*07:84",
                "HLA-B*07:85",
                "HLA-B*07:86",
                "HLA-B*07:87",
                "HLA-B*07:88",
                "HLA-B*07:89",
                "HLA-B*07:90",
                "HLA-B*07:91",
                "HLA-B*07:92",
                "HLA-B*07:93",
                "HLA-B*07:94",
                "HLA-B*07:95",
                "HLA-B*07:96",
                "HLA-B*07:97",
                "HLA-B*07:98",
                "HLA-B*07:99",
                "HLA-B*08:01",
                "HLA-B*08:02",
                "HLA-B*08:03",
                "HLA-B*08:04",
                "HLA-B*08:05",
                "HLA-B*08:07",
                "HLA-B*08:09",
                "HLA-B*08:10",
                "HLA-B*08:11",
                "HLA-B*08:12",
                "HLA-B*08:13",
                "HLA-B*08:14",
                "HLA-B*08:15",
                "HLA-B*08:16",
                "HLA-B*08:17",
                "HLA-B*08:18",
                "HLA-B*08:20",
                "HLA-B*08:21",
                "HLA-B*08:22",
                "HLA-B*08:23",
                "HLA-B*08:24",
                "HLA-B*08:25",
                "HLA-B*08:26",
                "HLA-B*08:27",
                "HLA-B*08:28",
                "HLA-B*08:29",
                "HLA-B*08:31",
                "HLA-B*08:32",
                "HLA-B*08:33",
                "HLA-B*08:34",
                "HLA-B*08:35",
                "HLA-B*08:36",
                "HLA-B*08:37",
                "HLA-B*08:38",
                "HLA-B*08:39",
                "HLA-B*08:40",
                "HLA-B*08:41",
                "HLA-B*08:42",
                "HLA-B*08:43",
                "HLA-B*08:44",
                "HLA-B*08:45",
                "HLA-B*08:46",
                "HLA-B*08:47",
                "HLA-B*08:48",
                "HLA-B*08:49",
                "HLA-B*08:50",
                "HLA-B*08:51",
                "HLA-B*08:52",
                "HLA-B*08:53",
                "HLA-B*08:54",
                "HLA-B*08:55",
                "HLA-B*08:56",
                "HLA-B*08:57",
                "HLA-B*08:58",
                "HLA-B*08:59",
                "HLA-B*08:60",
                "HLA-B*08:61",
                "HLA-B*08:62",
                "HLA-B*13:01",
                "HLA-B*13:02",
                "HLA-B*13:03",
                "HLA-B*13:04",
                "HLA-B*13:06",
                "HLA-B*13:09",
                "HLA-B*13:10",
                "HLA-B*13:11",
                "HLA-B*13:12",
                "HLA-B*13:13",
                "HLA-B*13:14",
                "HLA-B*13:15",
                "HLA-B*13:16",
                "HLA-B*13:17",
                "HLA-B*13:18",
                "HLA-B*13:19",
                "HLA-B*13:20",
                "HLA-B*13:21",
                "HLA-B*13:22",
                "HLA-B*13:23",
                "HLA-B*13:25",
                "HLA-B*13:26",
                "HLA-B*13:27",
                "HLA-B*13:28",
                "HLA-B*13:29",
                "HLA-B*13:30",
                "HLA-B*13:31",
                "HLA-B*13:32",
                "HLA-B*13:33",
                "HLA-B*13:34",
                "HLA-B*13:35",
                "HLA-B*13:36",
                "HLA-B*13:37",
                "HLA-B*13:38",
                "HLA-B*13:39",
                "HLA-B*14:01",
                "HLA-B*14:02",
                "HLA-B*14:03",
                "HLA-B*14:04",
                "HLA-B*14:05",
                "HLA-B*14:06",
                "HLA-B*14:08",
                "HLA-B*14:09",
                "HLA-B*14:10",
                "HLA-B*14:11",
                "HLA-B*14:12",
                "HLA-B*14:13",
                "HLA-B*14:14",
                "HLA-B*14:15",
                "HLA-B*14:16",
                "HLA-B*14:17",
                "HLA-B*14:18",
                "HLA-B*15:01",
                "HLA-B*15:02",
                "HLA-B*15:03",
                "HLA-B*15:04",
                "HLA-B*15:05",
                "HLA-B*15:06",
                "HLA-B*15:07",
                "HLA-B*15:08",
                "HLA-B*15:09",
                "HLA-B*15:10",
                "HLA-B*15:101",
                "HLA-B*15:102",
                "HLA-B*15:103",
                "HLA-B*15:104",
                "HLA-B*15:105",
                "HLA-B*15:106",
                "HLA-B*15:107",
                "HLA-B*15:108",
                "HLA-B*15:109",
                "HLA-B*15:11",
                "HLA-B*15:110",
                "HLA-B*15:112",
                "HLA-B*15:113",
                "HLA-B*15:114",
                "HLA-B*15:115",
                "HLA-B*15:116",
                "HLA-B*15:117",
                "HLA-B*15:118",
                "HLA-B*15:119",
                "HLA-B*15:12",
                "HLA-B*15:120",
                "HLA-B*15:121",
                "HLA-B*15:122",
                "HLA-B*15:123",
                "HLA-B*15:124",
                "HLA-B*15:125",
                "HLA-B*15:126",
                "HLA-B*15:127",
                "HLA-B*15:128",
                "HLA-B*15:129",
                "HLA-B*15:13",
                "HLA-B*15:131",
                "HLA-B*15:132",
                "HLA-B*15:133",
                "HLA-B*15:134",
                "HLA-B*15:135",
                "HLA-B*15:136",
                "HLA-B*15:137",
                "HLA-B*15:138",
                "HLA-B*15:139",
                "HLA-B*15:14",
                "HLA-B*15:140",
                "HLA-B*15:141",
                "HLA-B*15:142",
                "HLA-B*15:143",
                "HLA-B*15:144",
                "HLA-B*15:145",
                "HLA-B*15:146",
                "HLA-B*15:147",
                "HLA-B*15:148",
                "HLA-B*15:15",
                "HLA-B*15:150",
                "HLA-B*15:151",
                "HLA-B*15:152",
                "HLA-B*15:153",
                "HLA-B*15:154",
                "HLA-B*15:155",
                "HLA-B*15:156",
                "HLA-B*15:157",
                "HLA-B*15:158",
                "HLA-B*15:159",
                "HLA-B*15:16",
                "HLA-B*15:160",
                "HLA-B*15:161",
                "HLA-B*15:162",
                "HLA-B*15:163",
                "HLA-B*15:164",
                "HLA-B*15:165",
                "HLA-B*15:166",
                "HLA-B*15:167",
                "HLA-B*15:168",
                "HLA-B*15:169",
                "HLA-B*15:17",
                "HLA-B*15:170",
                "HLA-B*15:171",
                "HLA-B*15:172",
                "HLA-B*15:173",
                "HLA-B*15:174",
                "HLA-B*15:175",
                "HLA-B*15:176",
                "HLA-B*15:177",
                "HLA-B*15:178",
                "HLA-B*15:179",
                "HLA-B*15:18",
                "HLA-B*15:180",
                "HLA-B*15:183",
                "HLA-B*15:184",
                "HLA-B*15:185",
                "HLA-B*15:186",
                "HLA-B*15:187",
                "HLA-B*15:188",
                "HLA-B*15:189",
                "HLA-B*15:19",
                "HLA-B*15:191",
                "HLA-B*15:192",
                "HLA-B*15:193",
                "HLA-B*15:194",
                "HLA-B*15:195",
                "HLA-B*15:196",
                "HLA-B*15:197",
                "HLA-B*15:198",
                "HLA-B*15:199",
                "HLA-B*15:20",
                "HLA-B*15:200",
                "HLA-B*15:201",
                "HLA-B*15:202",
                "HLA-B*15:21",
                "HLA-B*15:23",
                "HLA-B*15:24",
                "HLA-B*15:25",
                "HLA-B*15:27",
                "HLA-B*15:28",
                "HLA-B*15:29",
                "HLA-B*15:30",
                "HLA-B*15:31",
                "HLA-B*15:32",
                "HLA-B*15:33",
                "HLA-B*15:34",
                "HLA-B*15:35",
                "HLA-B*15:36",
                "HLA-B*15:37",
                "HLA-B*15:38",
                "HLA-B*15:39",
                "HLA-B*15:40",
                "HLA-B*15:42",
                "HLA-B*15:43",
                "HLA-B*15:44",
                "HLA-B*15:45",
                "HLA-B*15:46",
                "HLA-B*15:47",
                "HLA-B*15:48",
                "HLA-B*15:49",
                "HLA-B*15:50",
                "HLA-B*15:51",
                "HLA-B*15:52",
                "HLA-B*15:53",
                "HLA-B*15:54",
                "HLA-B*15:55",
                "HLA-B*15:56",
                "HLA-B*15:57",
                "HLA-B*15:58",
                "HLA-B*15:60",
                "HLA-B*15:61",
                "HLA-B*15:62",
                "HLA-B*15:63",
                "HLA-B*15:64",
                "HLA-B*15:65",
                "HLA-B*15:66",
                "HLA-B*15:67",
                "HLA-B*15:68",
                "HLA-B*15:69",
                "HLA-B*15:70",
                "HLA-B*15:71",
                "HLA-B*15:72",
                "HLA-B*15:73",
                "HLA-B*15:74",
                "HLA-B*15:75",
                "HLA-B*15:76",
                "HLA-B*15:77",
                "HLA-B*15:78",
                "HLA-B*15:80",
                "HLA-B*15:81",
                "HLA-B*15:82",
                "HLA-B*15:83",
                "HLA-B*15:84",
                "HLA-B*15:85",
                "HLA-B*15:86",
                "HLA-B*15:87",
                "HLA-B*15:88",
                "HLA-B*15:89",
                "HLA-B*15:90",
                "HLA-B*15:91",
                "HLA-B*15:92",
                "HLA-B*15:93",
                "HLA-B*15:95",
                "HLA-B*15:96",
                "HLA-B*15:97",
                "HLA-B*15:98",
                "HLA-B*15:99",
                "HLA-B*18:01",
                "HLA-B*18:02",
                "HLA-B*18:03",
                "HLA-B*18:04",
                "HLA-B*18:05",
                "HLA-B*18:06",
                "HLA-B*18:07",
                "HLA-B*18:08",
                "HLA-B*18:09",
                "HLA-B*18:10",
                "HLA-B*18:11",
                "HLA-B*18:12",
                "HLA-B*18:13",
                "HLA-B*18:14",
                "HLA-B*18:15",
                "HLA-B*18:18",
                "HLA-B*18:19",
                "HLA-B*18:20",
                "HLA-B*18:21",
                "HLA-B*18:22",
                "HLA-B*18:24",
                "HLA-B*18:25",
                "HLA-B*18:26",
                "HLA-B*18:27",
                "HLA-B*18:28",
                "HLA-B*18:29",
                "HLA-B*18:30",
                "HLA-B*18:31",
                "HLA-B*18:32",
                "HLA-B*18:33",
                "HLA-B*18:34",
                "HLA-B*18:35",
                "HLA-B*18:36",
                "HLA-B*18:37",
                "HLA-B*18:38",
                "HLA-B*18:39",
                "HLA-B*18:40",
                "HLA-B*18:41",
                "HLA-B*18:42",
                "HLA-B*18:43",
                "HLA-B*18:44",
                "HLA-B*18:45",
                "HLA-B*18:46",
                "HLA-B*18:47",
                "HLA-B*18:48",
                "HLA-B*18:49",
                "HLA-B*18:50",
                "HLA-B*27:01",
                "HLA-B*27:02",
                "HLA-B*27:03",
                "HLA-B*27:04",
                "HLA-B*27:05",
                "HLA-B*27:06",
                "HLA-B*27:07",
                "HLA-B*27:08",
                "HLA-B*27:09",
                "HLA-B*27:10",
                "HLA-B*27:11",
                "HLA-B*27:12",
                "HLA-B*27:13",
                "HLA-B*27:14",
                "HLA-B*27:15",
                "HLA-B*27:16",
                "HLA-B*27:17",
                "HLA-B*27:18",
                "HLA-B*27:19",
                "HLA-B*27:20",
                "HLA-B*27:21",
                "HLA-B*27:23",
                "HLA-B*27:24",
                "HLA-B*27:25",
                "HLA-B*27:26",
                "HLA-B*27:27",
                "HLA-B*27:28",
                "HLA-B*27:29",
                "HLA-B*27:30",
                "HLA-B*27:31",
                "HLA-B*27:32",
                "HLA-B*27:33",
                "HLA-B*27:34",
                "HLA-B*27:35",
                "HLA-B*27:36",
                "HLA-B*27:37",
                "HLA-B*27:38",
                "HLA-B*27:39",
                "HLA-B*27:40",
                "HLA-B*27:41",
                "HLA-B*27:42",
                "HLA-B*27:43",
                "HLA-B*27:44",
                "HLA-B*27:45",
                "HLA-B*27:46",
                "HLA-B*27:47",
                "HLA-B*27:48",
                "HLA-B*27:49",
                "HLA-B*27:50",
                "HLA-B*27:51",
                "HLA-B*27:52",
                "HLA-B*27:53",
                "HLA-B*27:54",
                "HLA-B*27:55",
                "HLA-B*27:56",
                "HLA-B*27:57",
                "HLA-B*27:58",
                "HLA-B*27:60",
                "HLA-B*27:61",
                "HLA-B*27:62",
                "HLA-B*27:63",
                "HLA-B*27:67",
                "HLA-B*27:68",
                "HLA-B*27:69",
                "HLA-B*35:01",
                "HLA-B*35:02",
                "HLA-B*35:03",
                "HLA-B*35:04",
                "HLA-B*35:05",
                "HLA-B*35:06",
                "HLA-B*35:07",
                "HLA-B*35:08",
                "HLA-B*35:09",
                "HLA-B*35:10",
                "HLA-B*35:100",
                "HLA-B*35:101",
                "HLA-B*35:102",
                "HLA-B*35:103",
                "HLA-B*35:104",
                "HLA-B*35:105",
                "HLA-B*35:106",
                "HLA-B*35:107",
                "HLA-B*35:108",
                "HLA-B*35:109",
                "HLA-B*35:11",
                "HLA-B*35:110",
                "HLA-B*35:111",
                "HLA-B*35:112",
                "HLA-B*35:113",
                "HLA-B*35:114",
                "HLA-B*35:115",
                "HLA-B*35:116",
                "HLA-B*35:117",
                "HLA-B*35:118",
                "HLA-B*35:119",
                "HLA-B*35:12",
                "HLA-B*35:120",
                "HLA-B*35:121",
                "HLA-B*35:122",
                "HLA-B*35:123",
                "HLA-B*35:124",
                "HLA-B*35:125",
                "HLA-B*35:126",
                "HLA-B*35:127",
                "HLA-B*35:128",
                "HLA-B*35:13",
                "HLA-B*35:131",
                "HLA-B*35:132",
                "HLA-B*35:133",
                "HLA-B*35:135",
                "HLA-B*35:136",
                "HLA-B*35:137",
                "HLA-B*35:138",
                "HLA-B*35:139",
                "HLA-B*35:14",
                "HLA-B*35:140",
                "HLA-B*35:141",
                "HLA-B*35:142",
                "HLA-B*35:143",
                "HLA-B*35:144",
                "HLA-B*35:15",
                "HLA-B*35:16",
                "HLA-B*35:17",
                "HLA-B*35:18",
                "HLA-B*35:19",
                "HLA-B*35:20",
                "HLA-B*35:21",
                "HLA-B*35:22",
                "HLA-B*35:23",
                "HLA-B*35:24",
                "HLA-B*35:25",
                "HLA-B*35:26",
                "HLA-B*35:27",
                "HLA-B*35:28",
                "HLA-B*35:29",
                "HLA-B*35:30",
                "HLA-B*35:31",
                "HLA-B*35:32",
                "HLA-B*35:33",
                "HLA-B*35:34",
                "HLA-B*35:35",
                "HLA-B*35:36",
                "HLA-B*35:37",
                "HLA-B*35:38",
                "HLA-B*35:39",
                "HLA-B*35:41",
                "HLA-B*35:42",
                "HLA-B*35:43",
                "HLA-B*35:44",
                "HLA-B*35:45",
                "HLA-B*35:46",
                "HLA-B*35:47",
                "HLA-B*35:48",
                "HLA-B*35:49",
                "HLA-B*35:50",
                "HLA-B*35:51",
                "HLA-B*35:52",
                "HLA-B*35:54",
                "HLA-B*35:55",
                "HLA-B*35:56",
                "HLA-B*35:57",
                "HLA-B*35:58",
                "HLA-B*35:59",
                "HLA-B*35:60",
                "HLA-B*35:61",
                "HLA-B*35:62",
                "HLA-B*35:63",
                "HLA-B*35:64",
                "HLA-B*35:66",
                "HLA-B*35:67",
                "HLA-B*35:68",
                "HLA-B*35:69",
                "HLA-B*35:70",
                "HLA-B*35:71",
                "HLA-B*35:72",
                "HLA-B*35:74",
                "HLA-B*35:75",
                "HLA-B*35:76",
                "HLA-B*35:77",
                "HLA-B*35:78",
                "HLA-B*35:79",
                "HLA-B*35:80",
                "HLA-B*35:81",
                "HLA-B*35:82",
                "HLA-B*35:83",
                "HLA-B*35:84",
                "HLA-B*35:85",
                "HLA-B*35:86",
                "HLA-B*35:87",
                "HLA-B*35:88",
                "HLA-B*35:89",
                "HLA-B*35:90",
                "HLA-B*35:91",
                "HLA-B*35:92",
                "HLA-B*35:93",
                "HLA-B*35:94",
                "HLA-B*35:95",
                "HLA-B*35:96",
                "HLA-B*35:97",
                "HLA-B*35:98",
                "HLA-B*35:99",
                "HLA-B*37:01",
                "HLA-B*37:02",
                "HLA-B*37:04",
                "HLA-B*37:05",
                "HLA-B*37:06",
                "HLA-B*37:07",
                "HLA-B*37:08",
                "HLA-B*37:09",
                "HLA-B*37:10",
                "HLA-B*37:11",
                "HLA-B*37:12",
                "HLA-B*37:13",
                "HLA-B*37:14",
                "HLA-B*37:15",
                "HLA-B*37:17",
                "HLA-B*37:18",
                "HLA-B*37:19",
                "HLA-B*37:20",
                "HLA-B*37:21",
                "HLA-B*37:22",
                "HLA-B*37:23",
                "HLA-B*38:01",
                "HLA-B*38:02",
                "HLA-B*38:03",
                "HLA-B*38:04",
                "HLA-B*38:05",
                "HLA-B*38:06",
                "HLA-B*38:07",
                "HLA-B*38:08",
                "HLA-B*38:09",
                "HLA-B*38:10",
                "HLA-B*38:11",
                "HLA-B*38:12",
                "HLA-B*38:13",
                "HLA-B*38:14",
                "HLA-B*38:15",
                "HLA-B*38:16",
                "HLA-B*38:17",
                "HLA-B*38:18",
                "HLA-B*38:19",
                "HLA-B*38:20",
                "HLA-B*38:21",
                "HLA-B*38:22",
                "HLA-B*38:23",
                "HLA-B*39:01",
                "HLA-B*39:02",
                "HLA-B*39:03",
                "HLA-B*39:04",
                "HLA-B*39:05",
                "HLA-B*39:06",
                "HLA-B*39:07",
                "HLA-B*39:08",
                "HLA-B*39:09",
                "HLA-B*39:10",
                "HLA-B*39:11",
                "HLA-B*39:12",
                "HLA-B*39:13",
                "HLA-B*39:14",
                "HLA-B*39:15",
                "HLA-B*39:16",
                "HLA-B*39:17",
                "HLA-B*39:18",
                "HLA-B*39:19",
                "HLA-B*39:20",
                "HLA-B*39:22",
                "HLA-B*39:23",
                "HLA-B*39:24",
                "HLA-B*39:26",
                "HLA-B*39:27",
                "HLA-B*39:28",
                "HLA-B*39:29",
                "HLA-B*39:30",
                "HLA-B*39:31",
                "HLA-B*39:32",
                "HLA-B*39:33",
                "HLA-B*39:34",
                "HLA-B*39:35",
                "HLA-B*39:36",
                "HLA-B*39:37",
                "HLA-B*39:39",
                "HLA-B*39:41",
                "HLA-B*39:42",
                "HLA-B*39:43",
                "HLA-B*39:44",
                "HLA-B*39:45",
                "HLA-B*39:46",
                "HLA-B*39:47",
                "HLA-B*39:48",
                "HLA-B*39:49",
                "HLA-B*39:50",
                "HLA-B*39:51",
                "HLA-B*39:52",
                "HLA-B*39:53",
                "HLA-B*39:54",
                "HLA-B*39:55",
                "HLA-B*39:56",
                "HLA-B*39:57",
                "HLA-B*39:58",
                "HLA-B*39:59",
                "HLA-B*39:60",
                "HLA-B*40:01",
                "HLA-B*40:02",
                "HLA-B*40:03",
                "HLA-B*40:04",
                "HLA-B*40:05",
                "HLA-B*40:06",
                "HLA-B*40:07",
                "HLA-B*40:08",
                "HLA-B*40:09",
                "HLA-B*40:10",
                "HLA-B*40:100",
                "HLA-B*40:101",
                "HLA-B*40:102",
                "HLA-B*40:103",
                "HLA-B*40:104",
                "HLA-B*40:105",
                "HLA-B*40:106",
                "HLA-B*40:107",
                "HLA-B*40:108",
                "HLA-B*40:109",
                "HLA-B*40:11",
                "HLA-B*40:110",
                "HLA-B*40:111",
                "HLA-B*40:112",
                "HLA-B*40:113",
                "HLA-B*40:114",
                "HLA-B*40:115",
                "HLA-B*40:116",
                "HLA-B*40:117",
                "HLA-B*40:119",
                "HLA-B*40:12",
                "HLA-B*40:120",
                "HLA-B*40:121",
                "HLA-B*40:122",
                "HLA-B*40:123",
                "HLA-B*40:124",
                "HLA-B*40:125",
                "HLA-B*40:126",
                "HLA-B*40:127",
                "HLA-B*40:128",
                "HLA-B*40:129",
                "HLA-B*40:13",
                "HLA-B*40:130",
                "HLA-B*40:131",
                "HLA-B*40:132",
                "HLA-B*40:134",
                "HLA-B*40:135",
                "HLA-B*40:136",
                "HLA-B*40:137",
                "HLA-B*40:138",
                "HLA-B*40:139",
                "HLA-B*40:14",
                "HLA-B*40:140",
                "HLA-B*40:141",
                "HLA-B*40:143",
                "HLA-B*40:145",
                "HLA-B*40:146",
                "HLA-B*40:147",
                "HLA-B*40:15",
                "HLA-B*40:16",
                "HLA-B*40:18",
                "HLA-B*40:19",
                "HLA-B*40:20",
                "HLA-B*40:21",
                "HLA-B*40:23",
                "HLA-B*40:24",
                "HLA-B*40:25",
                "HLA-B*40:26",
                "HLA-B*40:27",
                "HLA-B*40:28",
                "HLA-B*40:29",
                "HLA-B*40:30",
                "HLA-B*40:31",
                "HLA-B*40:32",
                "HLA-B*40:33",
                "HLA-B*40:34",
                "HLA-B*40:35",
                "HLA-B*40:36",
                "HLA-B*40:37",
                "HLA-B*40:38",
                "HLA-B*40:39",
                "HLA-B*40:40",
                "HLA-B*40:42",
                "HLA-B*40:43",
                "HLA-B*40:44",
                "HLA-B*40:45",
                "HLA-B*40:46",
                "HLA-B*40:47",
                "HLA-B*40:48",
                "HLA-B*40:49",
                "HLA-B*40:50",
                "HLA-B*40:51",
                "HLA-B*40:52",
                "HLA-B*40:53",
                "HLA-B*40:54",
                "HLA-B*40:55",
                "HLA-B*40:56",
                "HLA-B*40:57",
                "HLA-B*40:58",
                "HLA-B*40:59",
                "HLA-B*40:60",
                "HLA-B*40:61",
                "HLA-B*40:62",
                "HLA-B*40:63",
                "HLA-B*40:64",
                "HLA-B*40:65",
                "HLA-B*40:66",
                "HLA-B*40:67",
                "HLA-B*40:68",
                "HLA-B*40:69",
                "HLA-B*40:70",
                "HLA-B*40:71",
                "HLA-B*40:72",
                "HLA-B*40:73",
                "HLA-B*40:74",
                "HLA-B*40:75",
                "HLA-B*40:76",
                "HLA-B*40:77",
                "HLA-B*40:78",
                "HLA-B*40:79",
                "HLA-B*40:80",
                "HLA-B*40:81",
                "HLA-B*40:82",
                "HLA-B*40:83",
                "HLA-B*40:84",
                "HLA-B*40:85",
                "HLA-B*40:86",
                "HLA-B*40:87",
                "HLA-B*40:88",
                "HLA-B*40:89",
                "HLA-B*40:90",
                "HLA-B*40:91",
                "HLA-B*40:92",
                "HLA-B*40:93",
                "HLA-B*40:94",
                "HLA-B*40:95",
                "HLA-B*40:96",
                "HLA-B*40:97",
                "HLA-B*40:98",
                "HLA-B*40:99",
                "HLA-B*41:01",
                "HLA-B*41:02",
                "HLA-B*41:03",
                "HLA-B*41:04",
                "HLA-B*41:05",
                "HLA-B*41:06",
                "HLA-B*41:07",
                "HLA-B*41:08",
                "HLA-B*41:09",
                "HLA-B*41:10",
                "HLA-B*41:11",
                "HLA-B*41:12",
                "HLA-B*42:01",
                "HLA-B*42:02",
                "HLA-B*42:04",
                "HLA-B*42:05",
                "HLA-B*42:06",
                "HLA-B*42:07",
                "HLA-B*42:08",
                "HLA-B*42:09",
                "HLA-B*42:10",
                "HLA-B*42:11",
                "HLA-B*42:12",
                "HLA-B*42:13",
                "HLA-B*42:14",
                "HLA-B*44:02",
                "HLA-B*44:03",
                "HLA-B*44:04",
                "HLA-B*44:05",
                "HLA-B*44:06",
                "HLA-B*44:07",
                "HLA-B*44:08",
                "HLA-B*44:09",
                "HLA-B*44:10",
                "HLA-B*44:100",
                "HLA-B*44:101",
                "HLA-B*44:102",
                "HLA-B*44:103",
                "HLA-B*44:104",
                "HLA-B*44:105",
                "HLA-B*44:106",
                "HLA-B*44:107",
                "HLA-B*44:109",
                "HLA-B*44:11",
                "HLA-B*44:110",
                "HLA-B*44:12",
                "HLA-B*44:13",
                "HLA-B*44:14",
                "HLA-B*44:15",
                "HLA-B*44:16",
                "HLA-B*44:17",
                "HLA-B*44:18",
                "HLA-B*44:20",
                "HLA-B*44:21",
                "HLA-B*44:22",
                "HLA-B*44:24",
                "HLA-B*44:25",
                "HLA-B*44:26",
                "HLA-B*44:27",
                "HLA-B*44:28",
                "HLA-B*44:29",
                "HLA-B*44:30",
                "HLA-B*44:31",
                "HLA-B*44:32",
                "HLA-B*44:33",
                "HLA-B*44:34",
                "HLA-B*44:35",
                "HLA-B*44:36",
                "HLA-B*44:37",
                "HLA-B*44:38",
                "HLA-B*44:39",
                "HLA-B*44:40",
                "HLA-B*44:41",
                "HLA-B*44:42",
                "HLA-B*44:43",
                "HLA-B*44:44",
                "HLA-B*44:45",
                "HLA-B*44:46",
                "HLA-B*44:47",
                "HLA-B*44:48",
                "HLA-B*44:49",
                "HLA-B*44:50",
                "HLA-B*44:51",
                "HLA-B*44:53",
                "HLA-B*44:54",
                "HLA-B*44:55",
                "HLA-B*44:57",
                "HLA-B*44:59",
                "HLA-B*44:60",
                "HLA-B*44:62",
                "HLA-B*44:63",
                "HLA-B*44:64",
                "HLA-B*44:65",
                "HLA-B*44:66",
                "HLA-B*44:67",
                "HLA-B*44:68",
                "HLA-B*44:69",
                "HLA-B*44:70",
                "HLA-B*44:71",
                "HLA-B*44:72",
                "HLA-B*44:73",
                "HLA-B*44:74",
                "HLA-B*44:75",
                "HLA-B*44:76",
                "HLA-B*44:77",
                "HLA-B*44:78",
                "HLA-B*44:79",
                "HLA-B*44:80",
                "HLA-B*44:81",
                "HLA-B*44:82",
                "HLA-B*44:83",
                "HLA-B*44:84",
                "HLA-B*44:85",
                "HLA-B*44:86",
                "HLA-B*44:87",
                "HLA-B*44:88",
                "HLA-B*44:89",
                "HLA-B*44:90",
                "HLA-B*44:91",
                "HLA-B*44:92",
                "HLA-B*44:93",
                "HLA-B*44:94",
                "HLA-B*44:95",
                "HLA-B*44:96",
                "HLA-B*44:97",
                "HLA-B*44:98",
                "HLA-B*44:99",
                "HLA-B*45:01",
                "HLA-B*45:02",
                "HLA-B*45:03",
                "HLA-B*45:04",
                "HLA-B*45:05",
                "HLA-B*45:06",
                "HLA-B*45:07",
                "HLA-B*45:08",
                "HLA-B*45:09",
                "HLA-B*45:10",
                "HLA-B*45:11",
                "HLA-B*45:12",
                "HLA-B*46:01",
                "HLA-B*46:02",
                "HLA-B*46:03",
                "HLA-B*46:04",
                "HLA-B*46:05",
                "HLA-B*46:06",
                "HLA-B*46:08",
                "HLA-B*46:09",
                "HLA-B*46:10",
                "HLA-B*46:11",
                "HLA-B*46:12",
                "HLA-B*46:13",
                "HLA-B*46:14",
                "HLA-B*46:16",
                "HLA-B*46:17",
                "HLA-B*46:18",
                "HLA-B*46:19",
                "HLA-B*46:20",
                "HLA-B*46:21",
                "HLA-B*46:22",
                "HLA-B*46:23",
                "HLA-B*46:24",
                "HLA-B*47:01",
                "HLA-B*47:02",
                "HLA-B*47:03",
                "HLA-B*47:04",
                "HLA-B*47:05",
                "HLA-B*47:06",
                "HLA-B*47:07",
                "HLA-B*48:01",
                "HLA-B*48:02",
                "HLA-B*48:03",
                "HLA-B*48:04",
                "HLA-B*48:05",
                "HLA-B*48:06",
                "HLA-B*48:07",
                "HLA-B*48:08",
                "HLA-B*48:09",
                "HLA-B*48:10",
                "HLA-B*48:11",
                "HLA-B*48:12",
                "HLA-B*48:13",
                "HLA-B*48:14",
                "HLA-B*48:15",
                "HLA-B*48:16",
                "HLA-B*48:17",
                "HLA-B*48:18",
                "HLA-B*48:19",
                "HLA-B*48:20",
                "HLA-B*48:21",
                "HLA-B*48:22",
                "HLA-B*48:23",
                "HLA-B*49:01",
                "HLA-B*49:02",
                "HLA-B*49:03",
                "HLA-B*49:04",
                "HLA-B*49:05",
                "HLA-B*49:06",
                "HLA-B*49:07",
                "HLA-B*49:08",
                "HLA-B*49:09",
                "HLA-B*49:10",
                "HLA-B*50:01",
                "HLA-B*50:02",
                "HLA-B*50:04",
                "HLA-B*50:05",
                "HLA-B*50:06",
                "HLA-B*50:07",
                "HLA-B*50:08",
                "HLA-B*50:09",
                "HLA-B*51:01",
                "HLA-B*51:02",
                "HLA-B*51:03",
                "HLA-B*51:04",
                "HLA-B*51:05",
                "HLA-B*51:06",
                "HLA-B*51:07",
                "HLA-B*51:08",
                "HLA-B*51:09",
                "HLA-B*51:12",
                "HLA-B*51:13",
                "HLA-B*51:14",
                "HLA-B*51:15",
                "HLA-B*51:16",
                "HLA-B*51:17",
                "HLA-B*51:18",
                "HLA-B*51:19",
                "HLA-B*51:20",
                "HLA-B*51:21",
                "HLA-B*51:22",
                "HLA-B*51:23",
                "HLA-B*51:24",
                "HLA-B*51:26",
                "HLA-B*51:28",
                "HLA-B*51:29",
                "HLA-B*51:30",
                "HLA-B*51:31",
                "HLA-B*51:32",
                "HLA-B*51:33",
                "HLA-B*51:34",
                "HLA-B*51:35",
                "HLA-B*51:36",
                "HLA-B*51:37",
                "HLA-B*51:38",
                "HLA-B*51:39",
                "HLA-B*51:40",
                "HLA-B*51:42",
                "HLA-B*51:43",
                "HLA-B*51:45",
                "HLA-B*51:46",
                "HLA-B*51:48",
                "HLA-B*51:49",
                "HLA-B*51:50",
                "HLA-B*51:51",
                "HLA-B*51:52",
                "HLA-B*51:53",
                "HLA-B*51:54",
                "HLA-B*51:55",
                "HLA-B*51:56",
                "HLA-B*51:57",
                "HLA-B*51:58",
                "HLA-B*51:59",
                "HLA-B*51:60",
                "HLA-B*51:61",
                "HLA-B*51:62",
                "HLA-B*51:63",
                "HLA-B*51:64",
                "HLA-B*51:65",
                "HLA-B*51:66",
                "HLA-B*51:67",
                "HLA-B*51:68",
                "HLA-B*51:69",
                "HLA-B*51:70",
                "HLA-B*51:71",
                "HLA-B*51:72",
                "HLA-B*51:73",
                "HLA-B*51:74",
                "HLA-B*51:75",
                "HLA-B*51:76",
                "HLA-B*51:77",
                "HLA-B*51:78",
                "HLA-B*51:79",
                "HLA-B*51:80",
                "HLA-B*51:81",
                "HLA-B*51:82",
                "HLA-B*51:83",
                "HLA-B*51:84",
                "HLA-B*51:85",
                "HLA-B*51:86",
                "HLA-B*51:87",
                "HLA-B*51:88",
                "HLA-B*51:89",
                "HLA-B*51:90",
                "HLA-B*51:91",
                "HLA-B*51:92",
                "HLA-B*51:93",
                "HLA-B*51:94",
                "HLA-B*51:95",
                "HLA-B*51:96",
                "HLA-B*52:01",
                "HLA-B*52:02",
                "HLA-B*52:03",
                "HLA-B*52:04",
                "HLA-B*52:05",
                "HLA-B*52:06",
                "HLA-B*52:07",
                "HLA-B*52:08",
                "HLA-B*52:09",
                "HLA-B*52:10",
                "HLA-B*52:11",
                "HLA-B*52:12",
                "HLA-B*52:13",
                "HLA-B*52:14",
                "HLA-B*52:15",
                "HLA-B*52:16",
                "HLA-B*52:17",
                "HLA-B*52:18",
                "HLA-B*52:19",
                "HLA-B*52:20",
                "HLA-B*52:21",
                "HLA-B*53:01",
                "HLA-B*53:02",
                "HLA-B*53:03",
                "HLA-B*53:04",
                "HLA-B*53:05",
                "HLA-B*53:06",
                "HLA-B*53:07",
                "HLA-B*53:08",
                "HLA-B*53:09",
                "HLA-B*53:10",
                "HLA-B*53:11",
                "HLA-B*53:12",
                "HLA-B*53:13",
                "HLA-B*53:14",
                "HLA-B*53:15",
                "HLA-B*53:16",
                "HLA-B*53:17",
                "HLA-B*53:18",
                "HLA-B*53:19",
                "HLA-B*53:20",
                "HLA-B*53:21",
                "HLA-B*53:22",
                "HLA-B*53:23",
                "HLA-B*54:01",
                "HLA-B*54:02",
                "HLA-B*54:03",
                "HLA-B*54:04",
                "HLA-B*54:06",
                "HLA-B*54:07",
                "HLA-B*54:09",
                "HLA-B*54:10",
                "HLA-B*54:11",
                "HLA-B*54:12",
                "HLA-B*54:13",
                "HLA-B*54:14",
                "HLA-B*54:15",
                "HLA-B*54:16",
                "HLA-B*54:17",
                "HLA-B*54:18",
                "HLA-B*54:19",
                "HLA-B*54:20",
                "HLA-B*54:21",
                "HLA-B*54:22",
                "HLA-B*54:23",
                "HLA-B*55:01",
                "HLA-B*55:02",
                "HLA-B*55:03",
                "HLA-B*55:04",
                "HLA-B*55:05",
                "HLA-B*55:07",
                "HLA-B*55:08",
                "HLA-B*55:09",
                "HLA-B*55:10",
                "HLA-B*55:11",
                "HLA-B*55:12",
                "HLA-B*55:13",
                "HLA-B*55:14",
                "HLA-B*55:15",
                "HLA-B*55:16",
                "HLA-B*55:17",
                "HLA-B*55:18",
                "HLA-B*55:19",
                "HLA-B*55:20",
                "HLA-B*55:21",
                "HLA-B*55:22",
                "HLA-B*55:23",
                "HLA-B*55:24",
                "HLA-B*55:25",
                "HLA-B*55:26",
                "HLA-B*55:27",
                "HLA-B*55:28",
                "HLA-B*55:29",
                "HLA-B*55:30",
                "HLA-B*55:31",
                "HLA-B*55:32",
                "HLA-B*55:33",
                "HLA-B*55:34",
                "HLA-B*55:35",
                "HLA-B*55:36",
                "HLA-B*55:37",
                "HLA-B*55:38",
                "HLA-B*55:39",
                "HLA-B*55:40",
                "HLA-B*55:41",
                "HLA-B*55:42",
                "HLA-B*55:43",
                "HLA-B*56:01",
                "HLA-B*56:02",
                "HLA-B*56:03",
                "HLA-B*56:04",
                "HLA-B*56:05",
                "HLA-B*56:06",
                "HLA-B*56:07",
                "HLA-B*56:08",
                "HLA-B*56:09",
                "HLA-B*56:10",
                "HLA-B*56:11",
                "HLA-B*56:12",
                "HLA-B*56:13",
                "HLA-B*56:14",
                "HLA-B*56:15",
                "HLA-B*56:16",
                "HLA-B*56:17",
                "HLA-B*56:18",
                "HLA-B*56:20",
                "HLA-B*56:21",
                "HLA-B*56:22",
                "HLA-B*56:23",
                "HLA-B*56:24",
                "HLA-B*56:25",
                "HLA-B*56:26",
                "HLA-B*56:27",
                "HLA-B*56:29",
                "HLA-B*57:01",
                "HLA-B*57:02",
                "HLA-B*57:03",
                "HLA-B*57:04",
                "HLA-B*57:05",
                "HLA-B*57:06",
                "HLA-B*57:07",
                "HLA-B*57:08",
                "HLA-B*57:09",
                "HLA-B*57:10",
                "HLA-B*57:11",
                "HLA-B*57:12",
                "HLA-B*57:13",
                "HLA-B*57:14",
                "HLA-B*57:15",
                "HLA-B*57:16",
                "HLA-B*57:17",
                "HLA-B*57:18",
                "HLA-B*57:19",
                "HLA-B*57:20",
                "HLA-B*57:21",
                "HLA-B*57:22",
                "HLA-B*57:23",
                "HLA-B*57:24",
                "HLA-B*57:25",
                "HLA-B*57:26",
                "HLA-B*57:27",
                "HLA-B*57:29",
                "HLA-B*57:30",
                "HLA-B*57:31",
                "HLA-B*57:32",
                "HLA-B*58:01",
                "HLA-B*58:02",
                "HLA-B*58:04",
                "HLA-B*58:05",
                "HLA-B*58:06",
                "HLA-B*58:07",
                "HLA-B*58:08",
                "HLA-B*58:09",
                "HLA-B*58:11",
                "HLA-B*58:12",
                "HLA-B*58:13",
                "HLA-B*58:14",
                "HLA-B*58:15",
                "HLA-B*58:16",
                "HLA-B*58:18",
                "HLA-B*58:19",
                "HLA-B*58:20",
                "HLA-B*58:21",
                "HLA-B*58:22",
                "HLA-B*58:23",
                "HLA-B*58:24",
                "HLA-B*58:25",
                "HLA-B*58:26",
                "HLA-B*58:27",
                "HLA-B*58:28",
                "HLA-B*58:29",
                "HLA-B*58:30",
                "HLA-B*59:01",
                "HLA-B*59:02",
                "HLA-B*59:03",
                "HLA-B*59:04",
                "HLA-B*59:05",
                "HLA-B*67:01",
                "HLA-B*67:02",
                "HLA-B*73:01",
                "HLA-B*73:02",
                "HLA-B*78:01",
                "HLA-B*78:02",
                "HLA-B*78:03",
                "HLA-B*78:04",
                "HLA-B*78:05",
                "HLA-B*78:06",
                "HLA-B*78:07",
                "HLA-B*81:01",
                "HLA-B*81:02",
                "HLA-B*81:03",
                "HLA-B*81:05",
                "HLA-B*82:01",
                "HLA-B*82:02",
                "HLA-B*82:03",
                "HLA-B*83:01",
                "HLA-C*01:02",
                "HLA-C*01:03",
                "HLA-C*01:04",
                "HLA-C*01:05",
                "HLA-C*01:06",
                "HLA-C*01:07",
                "HLA-C*01:08",
                "HLA-C*01:09",
                "HLA-C*01:10",
                "HLA-C*01:11",
                "HLA-C*01:12",
                "HLA-C*01:13",
                "HLA-C*01:14",
                "HLA-C*01:15",
                "HLA-C*01:16",
                "HLA-C*01:17",
                "HLA-C*01:18",
                "HLA-C*01:19",
                "HLA-C*01:20",
                "HLA-C*01:21",
                "HLA-C*01:22",
                "HLA-C*01:23",
                "HLA-C*01:24",
                "HLA-C*01:25",
                "HLA-C*01:26",
                "HLA-C*01:27",
                "HLA-C*01:28",
                "HLA-C*01:29",
                "HLA-C*01:30",
                "HLA-C*01:31",
                "HLA-C*01:32",
                "HLA-C*01:33",
                "HLA-C*01:34",
                "HLA-C*01:35",
                "HLA-C*01:36",
                "HLA-C*01:38",
                "HLA-C*01:39",
                "HLA-C*01:40",
                "HLA-C*02:02",
                "HLA-C*02:03",
                "HLA-C*02:04",
                "HLA-C*02:05",
                "HLA-C*02:06",
                "HLA-C*02:07",
                "HLA-C*02:08",
                "HLA-C*02:09",
                "HLA-C*02:10",
                "HLA-C*02:11",
                "HLA-C*02:12",
                "HLA-C*02:13",
                "HLA-C*02:14",
                "HLA-C*02:15",
                "HLA-C*02:16",
                "HLA-C*02:17",
                "HLA-C*02:18",
                "HLA-C*02:19",
                "HLA-C*02:20",
                "HLA-C*02:21",
                "HLA-C*02:22",
                "HLA-C*02:23",
                "HLA-C*02:24",
                "HLA-C*02:26",
                "HLA-C*02:27",
                "HLA-C*02:28",
                "HLA-C*02:29",
                "HLA-C*02:30",
                "HLA-C*02:31",
                "HLA-C*02:32",
                "HLA-C*02:33",
                "HLA-C*02:34",
                "HLA-C*02:35",
                "HLA-C*02:36",
                "HLA-C*02:37",
                "HLA-C*02:39",
                "HLA-C*02:40",
                "HLA-C*03:01",
                "HLA-C*03:02",
                "HLA-C*03:03",
                "HLA-C*03:04",
                "HLA-C*03:05",
                "HLA-C*03:06",
                "HLA-C*03:07",
                "HLA-C*03:08",
                "HLA-C*03:09",
                "HLA-C*03:10",
                "HLA-C*03:11",
                "HLA-C*03:12",
                "HLA-C*03:13",
                "HLA-C*03:14",
                "HLA-C*03:15",
                "HLA-C*03:16",
                "HLA-C*03:17",
                "HLA-C*03:18",
                "HLA-C*03:19",
                "HLA-C*03:21",
                "HLA-C*03:23",
                "HLA-C*03:24",
                "HLA-C*03:25",
                "HLA-C*03:26",
                "HLA-C*03:27",
                "HLA-C*03:28",
                "HLA-C*03:29",
                "HLA-C*03:30",
                "HLA-C*03:31",
                "HLA-C*03:32",
                "HLA-C*03:33",
                "HLA-C*03:34",
                "HLA-C*03:35",
                "HLA-C*03:36",
                "HLA-C*03:37",
                "HLA-C*03:38",
                "HLA-C*03:39",
                "HLA-C*03:40",
                "HLA-C*03:41",
                "HLA-C*03:42",
                "HLA-C*03:43",
                "HLA-C*03:44",
                "HLA-C*03:45",
                "HLA-C*03:46",
                "HLA-C*03:47",
                "HLA-C*03:48",
                "HLA-C*03:49",
                "HLA-C*03:50",
                "HLA-C*03:51",
                "HLA-C*03:52",
                "HLA-C*03:53",
                "HLA-C*03:54",
                "HLA-C*03:55",
                "HLA-C*03:56",
                "HLA-C*03:57",
                "HLA-C*03:58",
                "HLA-C*03:59",
                "HLA-C*03:60",
                "HLA-C*03:61",
                "HLA-C*03:62",
                "HLA-C*03:63",
                "HLA-C*03:64",
                "HLA-C*03:65",
                "HLA-C*03:66",
                "HLA-C*03:67",
                "HLA-C*03:68",
                "HLA-C*03:69",
                "HLA-C*03:70",
                "HLA-C*03:71",
                "HLA-C*03:72",
                "HLA-C*03:73",
                "HLA-C*03:74",
                "HLA-C*03:75",
                "HLA-C*03:76",
                "HLA-C*03:77",
                "HLA-C*03:78",
                "HLA-C*03:79",
                "HLA-C*03:80",
                "HLA-C*03:81",
                "HLA-C*03:82",
                "HLA-C*03:83",
                "HLA-C*03:84",
                "HLA-C*03:85",
                "HLA-C*03:86",
                "HLA-C*03:87",
                "HLA-C*03:88",
                "HLA-C*03:89",
                "HLA-C*03:90",
                "HLA-C*03:91",
                "HLA-C*03:92",
                "HLA-C*03:93",
                "HLA-C*03:94",
                "HLA-C*04:01",
                "HLA-C*04:03",
                "HLA-C*04:04",
                "HLA-C*04:05",
                "HLA-C*04:06",
                "HLA-C*04:07",
                "HLA-C*04:08",
                "HLA-C*04:10",
                "HLA-C*04:11",
                "HLA-C*04:12",
                "HLA-C*04:13",
                "HLA-C*04:14",
                "HLA-C*04:15",
                "HLA-C*04:16",
                "HLA-C*04:17",
                "HLA-C*04:18",
                "HLA-C*04:19",
                "HLA-C*04:20",
                "HLA-C*04:23",
                "HLA-C*04:24",
                "HLA-C*04:25",
                "HLA-C*04:26",
                "HLA-C*04:27",
                "HLA-C*04:28",
                "HLA-C*04:29",
                "HLA-C*04:30",
                "HLA-C*04:31",
                "HLA-C*04:32",
                "HLA-C*04:33",
                "HLA-C*04:34",
                "HLA-C*04:35",
                "HLA-C*04:36",
                "HLA-C*04:37",
                "HLA-C*04:38",
                "HLA-C*04:39",
                "HLA-C*04:40",
                "HLA-C*04:41",
                "HLA-C*04:42",
                "HLA-C*04:43",
                "HLA-C*04:44",
                "HLA-C*04:45",
                "HLA-C*04:46",
                "HLA-C*04:47",
                "HLA-C*04:48",
                "HLA-C*04:49",
                "HLA-C*04:50",
                "HLA-C*04:51",
                "HLA-C*04:52",
                "HLA-C*04:53",
                "HLA-C*04:54",
                "HLA-C*04:55",
                "HLA-C*04:56",
                "HLA-C*04:57",
                "HLA-C*04:58",
                "HLA-C*04:60",
                "HLA-C*04:61",
                "HLA-C*04:62",
                "HLA-C*04:63",
                "HLA-C*04:64",
                "HLA-C*04:65",
                "HLA-C*04:66",
                "HLA-C*04:67",
                "HLA-C*04:68",
                "HLA-C*04:69",
                "HLA-C*04:70",
                "HLA-C*05:01",
                "HLA-C*05:03",
                "HLA-C*05:04",
                "HLA-C*05:05",
                "HLA-C*05:06",
                "HLA-C*05:08",
                "HLA-C*05:09",
                "HLA-C*05:10",
                "HLA-C*05:11",
                "HLA-C*05:12",
                "HLA-C*05:13",
                "HLA-C*05:14",
                "HLA-C*05:15",
                "HLA-C*05:16",
                "HLA-C*05:17",
                "HLA-C*05:18",
                "HLA-C*05:19",
                "HLA-C*05:20",
                "HLA-C*05:21",
                "HLA-C*05:22",
                "HLA-C*05:23",
                "HLA-C*05:24",
                "HLA-C*05:25",
                "HLA-C*05:26",
                "HLA-C*05:27",
                "HLA-C*05:28",
                "HLA-C*05:29",
                "HLA-C*05:30",
                "HLA-C*05:31",
                "HLA-C*05:32",
                "HLA-C*05:33",
                "HLA-C*05:34",
                "HLA-C*05:35",
                "HLA-C*05:36",
                "HLA-C*05:37",
                "HLA-C*05:38",
                "HLA-C*05:39",
                "HLA-C*05:40",
                "HLA-C*05:41",
                "HLA-C*05:42",
                "HLA-C*05:43",
                "HLA-C*05:44",
                "HLA-C*05:45",
                "HLA-C*06:02",
                "HLA-C*06:03",
                "HLA-C*06:04",
                "HLA-C*06:05",
                "HLA-C*06:06",
                "HLA-C*06:07",
                "HLA-C*06:08",
                "HLA-C*06:09",
                "HLA-C*06:10",
                "HLA-C*06:11",
                "HLA-C*06:12",
                "HLA-C*06:13",
                "HLA-C*06:14",
                "HLA-C*06:15",
                "HLA-C*06:17",
                "HLA-C*06:18",
                "HLA-C*06:19",
                "HLA-C*06:20",
                "HLA-C*06:21",
                "HLA-C*06:22",
                "HLA-C*06:23",
                "HLA-C*06:24",
                "HLA-C*06:25",
                "HLA-C*06:26",
                "HLA-C*06:27",
                "HLA-C*06:28",
                "HLA-C*06:29",
                "HLA-C*06:30",
                "HLA-C*06:31",
                "HLA-C*06:32",
                "HLA-C*06:33",
                "HLA-C*06:34",
                "HLA-C*06:35",
                "HLA-C*06:36",
                "HLA-C*06:37",
                "HLA-C*06:38",
                "HLA-C*06:39",
                "HLA-C*06:40",
                "HLA-C*06:41",
                "HLA-C*06:42",
                "HLA-C*06:43",
                "HLA-C*06:44",
                "HLA-C*06:45",
                "HLA-C*07:01",
                "HLA-C*07:02",
                "HLA-C*07:03",
                "HLA-C*07:04",
                "HLA-C*07:05",
                "HLA-C*07:06",
                "HLA-C*07:07",
                "HLA-C*07:08",
                "HLA-C*07:09",
                "HLA-C*07:10",
                "HLA-C*07:100",
                "HLA-C*07:101",
                "HLA-C*07:102",
                "HLA-C*07:103",
                "HLA-C*07:105",
                "HLA-C*07:106",
                "HLA-C*07:107",
                "HLA-C*07:108",
                "HLA-C*07:109",
                "HLA-C*07:11",
                "HLA-C*07:110",
                "HLA-C*07:111",
                "HLA-C*07:112",
                "HLA-C*07:113",
                "HLA-C*07:114",
                "HLA-C*07:115",
                "HLA-C*07:116",
                "HLA-C*07:117",
                "HLA-C*07:118",
                "HLA-C*07:119",
                "HLA-C*07:12",
                "HLA-C*07:120",
                "HLA-C*07:122",
                "HLA-C*07:123",
                "HLA-C*07:124",
                "HLA-C*07:125",
                "HLA-C*07:126",
                "HLA-C*07:127",
                "HLA-C*07:128",
                "HLA-C*07:129",
                "HLA-C*07:13",
                "HLA-C*07:130",
                "HLA-C*07:131",
                "HLA-C*07:132",
                "HLA-C*07:133",
                "HLA-C*07:134",
                "HLA-C*07:135",
                "HLA-C*07:136",
                "HLA-C*07:137",
                "HLA-C*07:138",
                "HLA-C*07:139",
                "HLA-C*07:14",
                "HLA-C*07:140",
                "HLA-C*07:141",
                "HLA-C*07:142",
                "HLA-C*07:143",
                "HLA-C*07:144",
                "HLA-C*07:145",
                "HLA-C*07:146",
                "HLA-C*07:147",
                "HLA-C*07:148",
                "HLA-C*07:149",
                "HLA-C*07:15",
                "HLA-C*07:16",
                "HLA-C*07:17",
                "HLA-C*07:18",
                "HLA-C*07:19",
                "HLA-C*07:20",
                "HLA-C*07:21",
                "HLA-C*07:22",
                "HLA-C*07:23",
                "HLA-C*07:24",
                "HLA-C*07:25",
                "HLA-C*07:26",
                "HLA-C*07:27",
                "HLA-C*07:28",
                "HLA-C*07:29",
                "HLA-C*07:30",
                "HLA-C*07:31",
                "HLA-C*07:35",
                "HLA-C*07:36",
                "HLA-C*07:37",
                "HLA-C*07:38",
                "HLA-C*07:39",
                "HLA-C*07:40",
                "HLA-C*07:41",
                "HLA-C*07:42",
                "HLA-C*07:43",
                "HLA-C*07:44",
                "HLA-C*07:45",
                "HLA-C*07:46",
                "HLA-C*07:47",
                "HLA-C*07:48",
                "HLA-C*07:49",
                "HLA-C*07:50",
                "HLA-C*07:51",
                "HLA-C*07:52",
                "HLA-C*07:53",
                "HLA-C*07:54",
                "HLA-C*07:56",
                "HLA-C*07:57",
                "HLA-C*07:58",
                "HLA-C*07:59",
                "HLA-C*07:60",
                "HLA-C*07:62",
                "HLA-C*07:63",
                "HLA-C*07:64",
                "HLA-C*07:65",
                "HLA-C*07:66",
                "HLA-C*07:67",
                "HLA-C*07:68",
                "HLA-C*07:69",
                "HLA-C*07:70",
                "HLA-C*07:71",
                "HLA-C*07:72",
                "HLA-C*07:73",
                "HLA-C*07:74",
                "HLA-C*07:75",
                "HLA-C*07:76",
                "HLA-C*07:77",
                "HLA-C*07:78",
                "HLA-C*07:79",
                "HLA-C*07:80",
                "HLA-C*07:81",
                "HLA-C*07:82",
                "HLA-C*07:83",
                "HLA-C*07:84",
                "HLA-C*07:85",
                "HLA-C*07:86",
                "HLA-C*07:87",
                "HLA-C*07:88",
                "HLA-C*07:89",
                "HLA-C*07:90",
                "HLA-C*07:91",
                "HLA-C*07:92",
                "HLA-C*07:93",
                "HLA-C*07:94",
                "HLA-C*07:95",
                "HLA-C*07:96",
                "HLA-C*07:97",
                "HLA-C*07:99",
                "HLA-C*08:01",
                "HLA-C*08:02",
                "HLA-C*08:03",
                "HLA-C*08:04",
                "HLA-C*08:05",
                "HLA-C*08:06",
                "HLA-C*08:07",
                "HLA-C*08:08",
                "HLA-C*08:09",
                "HLA-C*08:10",
                "HLA-C*08:11",
                "HLA-C*08:12",
                "HLA-C*08:13",
                "HLA-C*08:14",
                "HLA-C*08:15",
                "HLA-C*08:16",
                "HLA-C*08:17",
                "HLA-C*08:18",
                "HLA-C*08:19",
                "HLA-C*08:20",
                "HLA-C*08:21",
                "HLA-C*08:22",
                "HLA-C*08:23",
                "HLA-C*08:24",
                "HLA-C*08:25",
                "HLA-C*08:27",
                "HLA-C*08:28",
                "HLA-C*08:29",
                "HLA-C*08:30",
                "HLA-C*08:31",
                "HLA-C*08:32",
                "HLA-C*08:33",
                "HLA-C*08:34",
                "HLA-C*08:35",
                "HLA-C*12:02",
                "HLA-C*12:03",
                "HLA-C*12:04",
                "HLA-C*12:05",
                "HLA-C*12:06",
                "HLA-C*12:07",
                "HLA-C*12:08",
                "HLA-C*12:09",
                "HLA-C*12:10",
                "HLA-C*12:11",
                "HLA-C*12:12",
                "HLA-C*12:13",
                "HLA-C*12:14",
                "HLA-C*12:15",
                "HLA-C*12:16",
                "HLA-C*12:17",
                "HLA-C*12:18",
                "HLA-C*12:19",
                "HLA-C*12:20",
                "HLA-C*12:21",
                "HLA-C*12:22",
                "HLA-C*12:23",
                "HLA-C*12:24",
                "HLA-C*12:25",
                "HLA-C*12:26",
                "HLA-C*12:27",
                "HLA-C*12:28",
                "HLA-C*12:29",
                "HLA-C*12:30",
                "HLA-C*12:31",
                "HLA-C*12:32",
                "HLA-C*12:33",
                "HLA-C*12:34",
                "HLA-C*12:35",
                "HLA-C*12:36",
                "HLA-C*12:37",
                "HLA-C*12:38",
                "HLA-C*12:40",
                "HLA-C*12:41",
                "HLA-C*12:43",
                "HLA-C*12:44",
                "HLA-C*14:02",
                "HLA-C*14:03",
                "HLA-C*14:04",
                "HLA-C*14:05",
                "HLA-C*14:06",
                "HLA-C*14:08",
                "HLA-C*14:09",
                "HLA-C*14:10",
                "HLA-C*14:11",
                "HLA-C*14:12",
                "HLA-C*14:13",
                "HLA-C*14:14",
                "HLA-C*14:15",
                "HLA-C*14:16",
                "HLA-C*14:17",
                "HLA-C*14:18",
                "HLA-C*14:19",
                "HLA-C*14:20",
                "HLA-C*15:02",
                "HLA-C*15:03",
                "HLA-C*15:04",
                "HLA-C*15:05",
                "HLA-C*15:06",
                "HLA-C*15:07",
                "HLA-C*15:08",
                "HLA-C*15:09",
                "HLA-C*15:10",
                "HLA-C*15:11",
                "HLA-C*15:12",
                "HLA-C*15:13",
                "HLA-C*15:15",
                "HLA-C*15:16",
                "HLA-C*15:17",
                "HLA-C*15:18",
                "HLA-C*15:19",
                "HLA-C*15:20",
                "HLA-C*15:21",
                "HLA-C*15:22",
                "HLA-C*15:23",
                "HLA-C*15:24",
                "HLA-C*15:25",
                "HLA-C*15:26",
                "HLA-C*15:27",
                "HLA-C*15:28",
                "HLA-C*15:29",
                "HLA-C*15:30",
                "HLA-C*15:31",
                "HLA-C*15:33",
                "HLA-C*15:34",
                "HLA-C*15:35",
                "HLA-C*16:01",
                "HLA-C*16:02",
                "HLA-C*16:04",
                "HLA-C*16:06",
                "HLA-C*16:07",
                "HLA-C*16:08",
                "HLA-C*16:09",
                "HLA-C*16:10",
                "HLA-C*16:11",
                "HLA-C*16:12",
                "HLA-C*16:13",
                "HLA-C*16:14",
                "HLA-C*16:15",
                "HLA-C*16:17",
                "HLA-C*16:18",
                "HLA-C*16:19",
                "HLA-C*16:20",
                "HLA-C*16:21",
                "HLA-C*16:22",
                "HLA-C*16:23",
                "HLA-C*16:24",
                "HLA-C*16:25",
                "HLA-C*16:26",
                "HLA-C*17:01",
                "HLA-C*17:02",
                "HLA-C*17:03",
                "HLA-C*17:04",
                "HLA-C*17:05",
                "HLA-C*17:06",
                "HLA-C*17:07",
                "HLA-C*18:01",
                "HLA-C*18:02",
                "HLA-C*18:03",
                "HLA-E*01:01",
                "HLA-G*01:01",
                "HLA-G*01:02",
                "HLA-G*01:03",
                "HLA-G*01:04",
                "HLA-G*01:06",
                "HLA-G*01:07",
                "HLA-G*01:08",
                "HLA-G*01:09",
            ]
        ),
    }

    if in_allele in alleles_dict[method]:
        return True
    else:
        return False


def combine_parsed_outputs(infile_list, outfile, top_score_metric, binding_cutoff=500):
    """"""
    fieldnames = []
    for input_file in infile_list:
        with open(input_file) as input_file_handle:
            reader = csv.DictReader(input_file_handle, delimiter="\t")
            if len(fieldnames) == 0:
                fieldnames = reader.fieldnames
            else:
                for fieldname in reader.fieldnames:
                    if fieldname not in fieldnames:
                        fieldnames.append(fieldname)

    rows = []
    for input_file in infile_list:
        with open(input_file) as input_file_handle:
            reader = csv.DictReader(input_file_handle, delimiter="\t")
            for row in reader:
                for fieldname in fieldnames:
                    if fieldname not in row:
                        row[fieldname] = "NA"
                rows.append(row)

    sorted_rows = sorted(rows, key=lambda row: (int(row["Sub-peptide Position"])))
    sorted_rows = sorted(
        sorted_rows,
        key=lambda row: (
            float(row["Corresponding Fold Change"])
            if row["Corresponding Fold Change"].isdigit()
            else float("inf")
        ),
        reverse=True,
    )

    if top_score_metric == "lowest":
        sorted_rows = sorted(
            sorted_rows,
            key=lambda row: (
                row["Gene Name"],
                row["Mutation"],
                float(row["Best MT Score"]),
            ),
        )
    elif top_score_metric == "median":
        sorted_rows = sorted(
            sorted_rows,
            key=lambda row: (
                row["Gene Name"],
                row["Mutation"],
                float(row["Median MT Score"]),
            ),
        )

    temp_dir = os.path.dirname(infile_list[0])
    tmp_file_name = os.path.join(temp_dir, f"{outfile}.tmp")
    tmp_file = open(tmp_file_name, "w")
    tsv_writer = csv.DictWriter(
        tmp_file, list(fieldnames), delimiter="\t", lineterminator="\n"
    )
    tsv_writer.writeheader()
    tsv_writer.writerows(sorted_rows)
    tmp_file.close()

    tmp_file = open(tmp_file_name)
    reader = csv.DictReader(tmp_file, delimiter="\t")
    fieldnames = reader.fieldnames
    outfile_obj = open(outfile, "w")
    writer = csv.DictWriter(
        outfile_obj, fieldnames, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()

    for entry in reader:
        if top_score_metric == "lowest":
            score = float(entry["Best MT Score"])
            fold_change = (
                sys.maxsize
                if entry["Corresponding Fold Change"] == "NA"
                else float(entry["Corresponding Fold Change"])
            )
        elif top_score_metric == "median":
            score = float(entry["Median MT Score"])
            fold_change = (
                sys.maxsize
                if entry["Median Fold Change"] == "NA"
                else float(entry["Median Fold Change"])
            )
        # if score > args.binding_threshold or fold_change < args.minimum_fold_change:
        if score > binding_cutoff:
            continue
        writer.writerow(entry)

    tmp_file.close()
    outfile_obj.close()


def add_ranking(infile, outfile, invcf, af_field="AB"):
    """Add the Ranking of Neoantigens based on binding affinity, fold change and allele frequency
    temporary solution, modify it later
    """

    indels_info = {}
    vcf_reader = vcf.Reader(open(f"{invcf}"))
    for record in vcf_reader:
        chrm = record.CHROM
        pos = record.POS
        ref = record.REF
        alt = str(record.ALT[0])
        af = float(record.INFO[af_field])
        # dp = int(record.INFO['DP'][0])
        indels_info[f"{chrm}:{pos}:{ref}:{alt}"] = {"AF": af}

    input_file_handle = open(infile)
    reader = csv.DictReader(input_file_handle, delimiter="\t")
    fieldnames = reader.fieldnames
    fieldnames.append("Ranking Score")

    rows = []
    for row in reader:
        index = "{}:{}:{}:{}".format(
            row["Chromosome"], row["Start"], row["Reference"], row["Variant"]
        )
        if index in indels_info:
            allele_freq = indels_info[index]["AF"]
        else:
            allele_freq = 0.0
        if row["Corresponding WT Score"] == "NA":
            wt_score = 100000.0
        else:
            wt_score = float(row["Corresponding WT Score"])

        row["Ranking Score"] = int(
            round(
                1.0 / float(row["Best MT Score"])
                + wt_score / float(row["Best MT Score"])
                + allele_freq * 100,
                2,
            )
        )
        rows.append(row)

    sorted_rows = sorted(
        rows, key=lambda row: (int(row["Ranking Score"])), reverse=True
    )

    tmpfile_name = f"{os.getpid()}.prerank.txt"
    tmp_file = open(tmpfile_name, "w")
    tsv_writer = csv.DictWriter(
        tmp_file, list(fieldnames), delimiter="\t", lineterminator="\n"
    )
    tsv_writer.writeheader()
    tsv_writer.writerows(sorted_rows)
    tmp_file.close()

    # Gene ranking
    gene_epitope_pairs = {}
    with open(tmpfile_name) as f:
        f.readline()
        for line in f:
            l = line.rstrip().split("\t")
            gene = l[10]
            epitope = l[15]
            neo_ranking = float(l[-1])
            gene_epitope_pairs[f"{gene}\t{epitope}"] = neo_ranking
    # key: gene, value: ranking scores
    gene_ranking = defaultdict(list)
    for i in gene_epitope_pairs:
        gene = i.split("\t")[0]
        neo_ranking = gene_epitope_pairs[i]
        gene_ranking[gene].append(neo_ranking)

    input_file_handle = open(tmpfile_name)
    reader = csv.DictReader(input_file_handle, delimiter="\t")
    fieldnames = reader.fieldnames
    fieldnames.append("Gene Ranking Score")
    rows = []
    for row in reader:
        row["Gene Ranking Score"] = int(np.median(gene_ranking[row["Gene Name"]]))
        rows.append(row)

    sorted_rows = sorted(
        rows, key=lambda row: (float(row["Ranking Score"])), reverse=True
    )

    out_file = open(outfile, "w")
    tsv_writer = csv.DictWriter(
        out_file, list(fieldnames), delimiter="\t", lineterminator="\n"
    )
    tsv_writer.writeheader()
    tsv_writer.writerows(sorted_rows)
    out_file.close()
    os.remove(tmpfile_name)


if __name__ == "__main__":
    result = OptiType_wrapper(sys.argv[1])
    alleles = optitype_parser(result)
    print(f"{os.path.basename(sys.argv[1])}\t{alleles}")
