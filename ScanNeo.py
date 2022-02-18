#!/usr/bin/env python
# ===============================================================================
__version__ = "2.2"
import sys
from pathlib import Path

root = str(Path(__file__).resolve().parents[0])
sys.path.append(root)
import re
import yaml
import os
import argparse
import numpy as np
import glob
from pyfaidx import Fasta
import subprocess
import secrets
import shutil
import datetime
import itertools
import vcf
from collections import OrderedDict, defaultdict
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
import tempfile
import multiprocessing
import textwrap
from io import StringIO

import ScanNeo_utils

MAX_WORKERS = multiprocessing.cpu_count()


def status_message(msg):
    print(msg)
    sys.stdout.flush()


chrms = {
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "chrX",
    "chrY",
}

chrms_dict = {
    "1": "chr1",
    "2": "chr2",
    "3": "chr3",
    "4": "chr4",
    "5": "chr5",
    "6": "chr6",
    "7": "chr7",
    "8": "chr8",
    "9": "chr9",
    "10": "chr10",
    "11": "chr11",
    "12": "chr12",
    "13": "chr13",
    "14": "chr14",
    "15": "chr15",
    "16": "chr16",
    "17": "chr17",
    "18": "chr18",
    "19": "chr19",
    "20": "chr20",
    "21": "chr21",
    "22": "chr22",
    "X": "chrX",
    "Y": "chrY",
    "MT": "chrM",
}


def external_tool_checking():
    """checking dependencies are installed"""
    software = [
        "picard",
        "vep",
        "sambamba",
        "bedtools",
        "bwa",
        "transIndel_build_RNA.py",
        "transIndel_call.py",
    ]
    cmd = "which"
    for each in software:
        try:
            path = subprocess.check_output([cmd, each], stderr=subprocess.STDOUT)
            path = str(path, "utf-8")
        except subprocess.CalledProcessError:
            print(
                "Checking for '" + each + "': ERROR - could not find '" + each + "'",
                file=sys.stderr,
            )
            print("Exiting.", file=sys.stderr)
            sys.exit(0)
        print("Checking for '" + each + "': found " + path)


def remove(infile):
    if os.path.isfile(infile):
        os.remove(infile)


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
        error_msg = f"Error happend!: {err}\n{err.output}"
    else:
        error_msg = ""
    if not error_msg:
        status_message(finish)
        return True
    else:
        status_message(error_msg)
        return False


config = ScanNeo_utils.config_getter()


def preprocessing(in_bam, ref="hg38", config=config):
    # discard splicing reads with alignment direction: XS:A:+ or XS:A:-
    name = os.path.splitext(in_bam)[0]
    if ref == "hg38":
        fasta = config["hg38_ref"]
    elif ref == "hg19":
        fasta = config["hg19_ref"]
    elif ref == "mm39":
        fasta = config["mm39_ref"]
    elif ref == "mm10":
        fasta = config["mm10_ref"]
    threads = config["threads"]

    # MAX_RECORDS_IN_RAM=500000
    duplication_remover = "picard MarkDuplicates I={0} O={1}.rmdup.bam M={1}.marked_dup_metrics.txt REMOVE_DUPLICATES=true CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT".format(
        in_bam, name
    )
    bam_filtering = 'sambamba view -f bam -t {1} -F "not (cigar =~ /N/ and [XS]!=null and not cigar =~ /I/ and not cigar =~ /D/) and mapping_quality > 0" {0}.rmdup.bam -o {0}.filter.bam'.format(
        name, threads
    )
    bam2fastq = "bedtools bamtofastq -i {0}.filter.bam -fq {0}.fq".format(name)

    bwa_mapping = "bwa mem -M -t {0} {1} {2}.fq | sambamba view -S -f bam /dev/stdin | sambamba sort -t {0} -o {2}.bwa.bam /dev/stdin".format(
        threads, fasta, name
    )
    index_bam = f"sambamba index {name}.bwa.bam"

    flag = False
    ret0 = run_cmd(duplication_remover, "Step0: duplicated reads removal!")
    if ret0:
        ret1 = run_cmd(bam_filtering, "Step1: splicing reads filtering finished!")
        if ret1:
            ret2 = run_cmd(bam2fastq, "Step2: FASTQ file generated!")
            if ret2:
                ret3 = run_cmd(bwa_mapping, "Step3: BWA-MEM mapping finished!")
                if ret3:
                    run_cmd(index_bam, "Step4: Indexing BWA BAM finished!")
                    flag = True

    remove(f"{name}.fq")
    # remove('{}.rmdup.bam'.format(name))
    # remove('{}.filter.bam'.format(name))
    if os.path.isfile(f"{name}.bwa.bam") and os.path.getsize(f"{name}.bwa.bam") > 0:
        return f"{name}.bwa.bam"
    else:
        sys.exit(f"Error: No {name}.bwa.bam generated!")
        return False


def scansv_caller(in_bam, ref="hg38", config=config):
    name = os.path.splitext(in_bam)[0].replace(".bwa", "")
    if ref == "hg38":
        fasta = config["hg38_ref"]
        annotation = config["hg38_anno"]
    elif ref == "hg19":
        fasta = config["hg19_ref"]
        annotation = config["hg19_anno"]
    elif ref == "mm39":
        fasta = config["mm39_ref"]
        annotation = config["hg39_anno"]
    elif ref == "mm10":
        fasta = config["mm10_ref"]
        annotation = config["mm10_anno"]

    scansv_build = f"transIndel_build_RNA.py -r {fasta} -g {annotation} -i {in_bam} -o {name}.indel.bam"
    scansv_call = "transIndel_call.py -i {0}.indel.bam -o {0} -d 10 -m 50 -l 1".format(
        name
    )

    flag = False
    ret1 = run_cmd(scansv_build, "Begin building indel.bam using transIndel,finished!")
    if ret1:
        ret2 = run_cmd(scansv_call, "Begin calling INDELs using transIndel,finished!")
        if ret2:
            flag = True
    if flag:
        remove(in_bam)
        remove(f"{in_bam}.bai")
        if (
            os.path.isfile(f"{name}.indel.vcf")
            and os.path.getsize(f"{name}.indel.vcf") > 0
        ):
            return f"{name}.indel.vcf"
        else:
            sys.exit(f"Error: No {name}.indel.vcf generated!")
            return False


# Function to reverse a string
def reverse(string):
    string = string[::-1]
    return string


def is_slippage(chrm, pos, ref, alt, indel_type):
    if indel_type == "INS":
        if repeat_checker(alt[1:]):
            pat = re.compile(rf"{alt[0]}{alt[1] * 4}")
        else:
            pat = re.compile(rf"{alt[0]}{alt[1:] * 4}")
        match = pat.match(genome_seq[chrm][pos - 1 : pos + 10])
        if match:
            return True
        else:
            return False
    elif indel_type == "DEL":
        if repeat_checker(ref[1:]):
            pat = re.compile(rf"{ref[0]}{ref[1] * 4}")
        else:
            pat = re.compile(rf"{ref[0]}{ref[1:] * 4}")
        match = pat.match(genome_seq[chrm][pos - 1 : pos + 10])
        if match:
            return True
        else:
            return False


def repeat_checker(s):
    if len(s) == s.count(s[0]):
        return True
    else:
        return False


def vcf_renewer(in_vcf, out_vcf, ref="hg38", slippage=False, config=config):
    """add REF and ALT sequences"""
    status_message(f"Begin deal with VCF: {in_vcf}")
    if ref == "hg38":
        fasta = config["hg38_ref"]
    elif ref == "hg19":
        fasta = config["hg19_ref"]
    elif ref == "mm39":
        fasta = config["mm39_ref"]
    elif ref == "mm10":
        fasta = config["mm10_ref"]

    genome_seq = Fasta(fasta, sequence_always_upper=True, as_raw=True)

    vcf_reader = vcf.Reader(open(f"{in_vcf}"))
    vcf_writer = vcf.Writer(open(f"{out_vcf}", "w"), vcf_reader)
    for record in vcf_reader:
        sv_type = record.INFO["SVTYPE"]
        chrm = record.CHROM
        pos = record.POS
        end = record.INFO["END"]
        if chrm in chrms:
            if sv_type == "INS":
                alt = str(record.ALT[0])
                if repeat_checker(alt[1:]):
                    pat = re.compile(rf"{alt[0]}{alt[1] * 4}")
                else:
                    pat = re.compile(rf"{alt[0]}{alt[1:] * 4}")
                match = pat.match(genome_seq[chrm][pos - 1 : pos + 10])
                if match:
                    is_slippage = True
                else:
                    is_slippage = False
            elif sv_type == "DEL":
                record.ALT = [genome_seq[chrm][pos - 1 : pos]]
                record.REF = genome_seq[chrm][pos - 1 : end]
                ref = record.REF
                if repeat_checker(ref[1:]):
                    pat = re.compile(rf"{ref[0]}{ref[1] * 4}")
                else:
                    pat = re.compile(rf"{ref[0]}{ref[1:] * 4}")
                match = pat.match(genome_seq[chrm][pos - 1 : pos + 10])
                if match:
                    is_slippage = True
                else:
                    is_slippage = False

            if not slippage:
                if not is_slippage:
                    vcf_writer.write_record(record)
            else:
                vcf_writer.write_record(record)

    vcf_writer.close()
    status_message("VCF add sequences [ref, alt] accomplished!")


def vep_caller(in_vcf, out_vcf, cutoff=0.01, ref="hg38"):
    if ref == "hg19":
        assembly = "GRCh37"
    elif ref == "hg38":
        assembly = "GRCh38"
    elif ref == "mm39":
        assembly = "GRCm39"
    elif ref == "mm10":
        assembly = "GRCm38"

    rnd_id = secrets.randbits(42)
    tmp_vcf = f"tmp.{rnd_id}.vcf"
    if ref.startswith("hg"):
        _af_gnomad = "--af_gnomad"
        _species = "homo_sapiens"
    else:
        _af_gnomad = ""
        _species = "mus_musculus"

    vep_cmd = f"vep --offline --force_overwrite --species {_species} --assembly {assembly} --input_file {in_vcf} --format vcf --output_file {tmp_vcf} \
        --vcf --symbol --terms SO --af {_af_gnomad} --plugin Downstream --plugin Wildtype --no_stats"

    ret = run_cmd(vep_cmd, "VEP begin annotation, VEP annotation finished!")
    flag = False
    if ret:
        status_message("VEP accomplished!")
        flag = True
    if flag:
        vcf_reader = vcf.Reader(open(f"{tmp_vcf}"))
        csq_format = ScanNeo_utils.parse_csq_format(vcf_reader)
        vcf_writer = vcf.Writer(open(f"{out_vcf}", "w"), vcf_reader)
        for record in vcf_reader:
            alleles_dict = ScanNeo_utils.resolve_alleles(record)
            alt = record.ALT[0]
            csq_allele = alleles_dict[alt]
            transcripts = ScanNeo_utils.parse_csq_entries_for_allele(
                record.INFO["CSQ"], csq_format, csq_allele
            )
            transcript_one = transcripts[0]
            AF = float(transcript_one["AF"]) if transcript_one["AF"] else 0.0
            try:
                gnomAD_AF = float(transcript_one["gnomAD_AF"])
            except KeyError:
                gnomAD_AF = 0.0
            except ValueError:
                gnomAD_AF = 0.0

            # vcf filtering allele frequency > 1% in 1000genomes or gnomAD database
            if AF > cutoff or gnomAD_AF > cutoff:
                pass
            else:
                vcf_writer.write_record(record)
        vcf_writer.close()
        status_message("VCF filtering accomplished!")
        remove(tmp_vcf)
        return out_vcf


################################################################################


def convert_vcf(input_file, temp_dir):
    print("Converting VCF to TSV")
    tsv_file = os.path.join(temp_dir, "tmp.tsv")
    convert_params = [
        input_file,
        tsv_file,
    ]
    ScanNeo_utils.convert_vcf(convert_params)
    print("Completed")


def generate_fasta(
    peptide_sequence_length, epitope_lengths, downstream_sequence_length, temp_dir
):
    print("Generating Variant Peptide FASTA and Key File")
    tsv_file = os.path.join(temp_dir, "tmp.tsv")
    fasta_file = os.path.join(temp_dir, "tmp.fasta")
    fasta_key_file = os.path.join(temp_dir, "tmp.fasta.key")
    generate_fasta_params = [
        tsv_file,
        str(peptide_sequence_length),
        str(max(epitope_lengths)),
        fasta_file,
        fasta_key_file,
        downstream_sequence_length,
    ]
    ScanNeo_utils.generate_fasta(generate_fasta_params)
    if os.path.isfile(fasta_file) and os.path.getsize(fasta_file) > 0:
        return True
        print("Completed")
    else:
        sys.exit("Error: No protein changing indels found in the VCF")


def generate_protein_fasta(
    input_vcf, peptide_sequence_length, epitope_lengths, downstream_length=1000
):
    if downstream_length == "full":
        downstream_sequence_length = 0
    elif str(downstream_length).isdigit():
        downstream_sequence_length = str(downstream_length)
    else:
        sys.exit(
            "The downstream sequence length needs to be a positive integer or 'full'"
        )

    temp_dir = tempfile.mkdtemp()
    convert_vcf(input_vcf, temp_dir)
    generate_fasta(
        peptide_sequence_length, epitope_lengths, downstream_sequence_length, temp_dir
    )
    return temp_dir


def parse_files(output_file, temp_dir):
    fasta_file = os.path.join(temp_dir, "tmp.fasta")
    fasta_key_file = os.path.join(temp_dir, "tmp.fasta.key")

    with open(fasta_key_file) as _fasta_key_file:
        keys = yaml.safe_load(_fasta_key_file)

    dataframe = OrderedDict()
    with open(fasta_file) as _fasta_file:
        for line in _fasta_file:
            key = line.rstrip().replace(">", "")
            sequence = _fasta_file.readline().rstrip()
            ids = keys[int(key)]
            for id in ids:
                (type, index) = id.split(".", 1)
                if index not in dataframe:
                    dataframe[index] = {}
                dataframe[index][type] = sequence

    with open(output_file, "w") as parsed_fasta_file:
        for index, sequences in dataframe.items():
            for type in ("WT", "MT"):
                parsed_fasta_file.write(f">{type}.{index}\n")
                parsed_fasta_file.write(f"{sequences[type]}\n")
    print("Completed")


def generate_protein_fasta_main(
    input_vcf,
    peptide_sequence_length,
    output_fasta_file,
    downstream_sequence_length=1000,
):
    temp_dir = generate_protein_fasta(
        input_vcf, peptide_sequence_length, "0", downstream_sequence_length
    )
    parse_files(output_fasta_file, temp_dir)


def iedb_caller(path_to_iedb, method, allele, epitope_length, temp_dir):
    """IEDB MHC caller (HLA class I)
    apply one method, one allele and one type of epitope length
    output: full path to iedb output file
    """
    methods = {"ann": "NetMHC", "netmhcpan": "NetMHCpan"}
    fasta = os.path.join(temp_dir, "tmp.fasta")

    if not allele_validator(allele, method):
        status_message(f"{allele} is not available for {methods[method]}")
        status_message("Skipping ...")
        return None

    status_message(
        f"Running IEDB on Allele {allele} and Epitope Length {epitope_length} with Method {methods[method]}"
    )

    iedb_mhc_i_executable = os.path.join(
        path_to_iedb, "mhc_i", "src", "predict_binding.py"
    )
    if not os.path.exists(iedb_mhc_i_executable):
        sys.exit("IEDB MHC I executable path doesn't exist %s" % iedb_mhc_i_executable)
    if not int(epitope_length) in {8, 9, 10, 11}:
        sys.exit("Epitope length should chose from {8,9,10,11}")

    iedb_out = os.path.join(
        temp_dir, ".".join([method, allele, str(epitope_length), "tsv"])
    )
    out_file = open(iedb_out, "w")

    cmd = f"{iedb_mhc_i_executable} {method} {allele} {epitope_length} {fasta}"
    response = subprocess.check_output(cmd, stdin=subprocess.PIPE, shell=True)
    # iedb_result = StringIO(response)
    # out_file.write(iedb_result.readline())
    # for line in iedb_result:
    #    tmp_l = line.rstrip('\n').split('\t')
    #    ic50  = float(tmp_l[-2])
    #    if ic50 <= binding_cutoff:
    #        out_file.write(line)
    out_file.write(response.decode("utf8"))
    out_file.close()
    status_message(
        f"Complete running IEDB on Allele {allele} and Epitope Length {epitope_length} with Method {methods[method]}"
    )
    return iedb_out


def call_iedb_and_parse_outputs(
    path_to_iedb,
    alleles,
    epitope_lengths,
    temp_dir,
    methods,
    top_score_metric,
    number_of_threads,
):
    """parallelized call iedb and combine allele epitope_length combination (per sample) to parsed files"""
    fasta_key_file = os.path.join(temp_dir, "tmp.fasta.key")
    tsv_file = os.path.join(temp_dir, "tmp.tsv")
    iterable = [methods, alleles, epitope_lengths]
    parameters = []
    for method, allele, epitope_length in itertools.product(*iterable):
        parameters.append((method, allele, epitope_length))

    output_files = []
    with ProcessPoolExecutor(max_workers=number_of_threads) as executor:
        to_do = {}
        for method, allele, epitope_length in parameters:
            # iedb_caller(path_to_iedb, method, allele, epitope_length, temp_dir):
            future = executor.submit(
                iedb_caller, path_to_iedb, method, allele, epitope_length, temp_dir
            )
            to_do[future] = f"{method}|{allele}|{epitope_length}"
        for future in as_completed(to_do):
            try:
                iedb_out = future.result()
                if iedb_out:
                    if os.path.isfile(iedb_out) and os.path.getsize(iedb_out) > 0:
                        output_files.append(iedb_out)
            except Exception as exc:
                error_msg = "error: " + str(exc)
                status_message(error_msg)
            else:
                error_msg = ""

    parsed_output_files = []
    if len(output_files) > 0:
        # store iedb output files with the same allele and epitope in a list
        allele_epitope_dict = defaultdict(list)
        for output_file in output_files:
            method, allele, epl, _ = os.path.basename(output_file).split(".")
            allele_epitope_dict[f"{allele}|{epl}"].append(output_file)

        with ProcessPoolExecutor(max_workers=number_of_threads) as executor:
            parse_to_do = {}

            for allele_epl in allele_epitope_dict:
                allele, epl = allele_epl.split("|")
                status_message(
                    f"Parsing IEDB Output for Allele {allele} and Epitope Length {epl}"
                )
                iedb_output_files = allele_epitope_dict[allele_epl]
                parsed_file = os.path.join(
                    temp_dir, ".".join([allele, str(epl), "parsed", "tsv"])
                )

                params = []
                params.extend(iedb_output_files)
                params.extend(
                    [tsv_file, fasta_key_file, parsed_file, "-m", top_score_metric]
                )

                future = executor.submit(ScanNeo_utils.parse_output, params)
                parse_to_do[future] = f"{allele}|{epl}"

            for future in as_completed(parse_to_do):
                try:
                    parsed_out = future.result()
                    allele, epl = parse_to_do[future].split("|")
                    if parsed_out:
                        if (
                            os.path.isfile(parsed_out)
                            and os.path.getsize(parsed_out) > 0
                        ):
                            parsed_output_files.append(parsed_out)
                except Exception as exc:
                    error_msg = "error: " + str(exc)
                    status_message(error_msg)
                else:
                    error_msg = ""
                    status_message(
                        f"Complete parsing IEDB Output for Allele {allele} and Epitope Length {epl}"
                    )
            status_message("Completed")

    return parsed_output_files


def combine_results(infile_list, outfile, top_score_metric, binding_cutoff):
    """combine prediction results and filter based on ic50 < 500nM"""
    if len(infile_list) > 0:
        ScanNeo_utils.combine_parsed_outputs(
            infile_list, outfile, top_score_metric, binding_cutoff
        )
        return outfile
    else:
        status_message("No valid allele available")
        status_message("No prediction results from IEDB")
        return None


def allele_validator(allele, method):
    """Input: an allele list
    Ouput: a valid allele frozenset
    """
    valid_allele_or_not = ScanNeo_utils.allele_checker(allele, method)
    return valid_allele_or_not


def OptiType_runner(in_bam):
    result = ScanNeo_utils.OptiType_wrapper(in_bam)
    alleles = ScanNeo_utils.optitype_parser(result)
    return alleles


def reference_proteome_filter(in_file):
    """TODO"""
    pass


def parse_args():
    parser = argparse.ArgumentParser(
        description="ScanNeo pipeline: neoantigen identification using RNA-seq data",
        epilog=textwrap.dedent(
            """Author: Ting-You Wang <tywang@umn.edu>, Hormel Institute, University of Minnesota, 2019"""
        ),
    )
    parser.add_argument(
        "-v", "--version", action="version", version=f"%(prog)s {__version__}"
    )
    sub_parsers = parser.add_subparsers(help="sub-command help", dest="sub_command")

    indel_parser = sub_parsers.add_parser(
        "indel",
        help="INDELs calling using transIndel",
        description="%(prog)s -i rnaseq_bam -r hg38",
        epilog=textwrap.dedent(
            """Author: Ting-You Wang <tywang@umn.edu>, Hormel Institute, University of Minnesota, 2019"""
        ),
    )
    indel_parser.add_argument(
        "-i",
        "--input",
        action="store",
        dest="input",
        help="RNA-seq alignment file (BAM)",
        required=True,
    )
    indel_parser.add_argument(
        "-r",
        "--ref",
        action="store",
        dest="ref",
        help="reference genome (default: %(default)s)",
        choices=["hg19", "hg38", "mm39", "mm10"],
        default="hg38",
    )

    vcf_parser = sub_parsers.add_parser(
        "anno",
        help="annotate INDELs (VCF) using VEP and filter based on 1KG and gnomAD",
        description="%(prog)s -i indel.vcf -c maf_cutoff -o output.vcf",
        epilog=textwrap.dedent(
            """Author: Ting-You Wang <tywang@umn.edu>, Hormel Institute, University of Minnesota, 2019"""
        ),
    )
    vcf_parser.add_argument(
        "-i",
        "--input",
        action="store",
        dest="input",
        help="input VCF file",
        required=True,
    )
    vcf_parser.add_argument(
        "-c",
        "--cutoff",
        action="store",
        dest="cutoff",
        help="MAF cutoff default: %(default)s",
        type=float,
        default=0.01,
    )
    vcf_parser.add_argument(
        "-s", "--slippage", action="store_true", dest="slippage", help="Keep slippage"
    )
    vcf_parser.add_argument(
        "-r",
        "--ref",
        action="store",
        dest="ref",
        help="reference genome (default: %(default)s)",
        choices=["hg19", "hg38", "mm39", "mm10"],
        default="hg38",
    )
    vcf_parser.add_argument(
        "-o",
        "--output",
        action="store",
        dest="output",
        help="output annotated and filtered vcf file (default: %(default)s)",
        default="output.vcf",
    )

    class ThreadAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            if values > MAX_WORKERS:
                parser.error(f"Maximum number of threads is {MAX_WORKERS}")
            setattr(namespace, self.dest, values)

    hla_parser = sub_parsers.add_parser(
        "hla",
        help="HLA genotyping and binding affinity prediction",
        description="%(prog)s -i vep.vcf --alleles allele1,allele2 -l peptide_sequence_length -o output.tsv",
        epilog=textwrap.dedent(
            """Author: Ting-You Wang <tywang@umn.edu>, Hormel Institute, University of Minnesota, 2019"""
        ),
    )
    hla_parser.add_argument(
        "-i",
        "--input",
        action="store",
        dest="vcf",
        help="VEP annotated and filtered VCF",
        required=True,
    )
    hla_parser.add_argument(
        "--alleles",
        action="store",
        dest="alleles",
        help="Name of the allele to use for epitope prediction. Multiple alleles can be specified using a comma-separated listinput HLA class I alleles",
    )
    hla_parser.add_argument(
        "-b",
        "--bam",
        action="store",
        dest="bam",
        help="Input RNA-Seq BAM file if you don't know sample HLA class I alleles. Works for human only!",
    )
    hla_parser.add_argument(
        "-l",
        "--length",
        action="store",
        dest="length",
        type=int,
        default=21,
        help="Length of the peptide sequence to use when creating the FASTA (default: %(default)s)",
    )
    hla_parser.add_argument(
        "-t",
        "--thread",
        action=ThreadAction,
        dest="thread",
        type=int,
        default=1,
        help="Number of threads (default: %(default)s)",
    )
    hla_parser.add_argument(
        "--af",
        action="store",
        dest="af_field",
        default="AB",
        help="The field name for allele frequency in VCF (default: %(default)s)",
    )
    hla_parser.add_argument(
        "--binding",
        action="store",
        dest="binding_threshold",
        type=int,
        default=500,
        help="binding threshold ic50 (default: %(default)s nM)",
    )
    hla_parser.add_argument(
        "-e",
        "--epitope-length",
        action="store",
        dest="epitope_lengths",
        help="Length of subpeptides (neoepitopes) to predict. Multiple epitope lengths can be specified using a comma-separated list. Typical epitope lengths vary between 8-11. (default: %(default)s)",
        default="8,9,10,11",
    )
    hla_parser.add_argument(
        "-p",
        "--path-to-iedb",
        action="store",
        dest="path_to_iedb",
        help="Directory that contains the local installation of IEDB",
        required=True,
    )
    hla_parser.add_argument(
        "-m",
        "--metric",
        action="store",
        dest="metric",
        choices=["lowest", "median"],
        default="lowest",
        help="The ic50 scoring metric to use when filtering epitopes by binding-threshold lowest: Best MT Score - lowest MT ic50 binding score of all chosen prediction methods. median: Median MT Score - median MT ic50 binding score of all chosen prediction methods. (default: %(default)s)",
    )
    hla_parser.add_argument(
        "-o",
        "--output",
        action="store",
        dest="output",
        help="output text file name, name.tsv",
    )

    fasta_parser = sub_parsers.add_parser(
        "fasta",
        help="Generate protein fasta from VEP-annotated VCF",
        description="%(prog)s -i vep.vcf -l peptide_sequence_length -o output.fasta",
        epilog=textwrap.dedent(
            """Author: Ting-You Wang <tywang@umn.edu>, Hormel Institute, University of Minnesota, 2019"""
        ),
    )
    fasta_parser.add_argument(
        "-i",
        "--input",
        action="store",
        dest="vcf",
        help="A VEP-annotated single-sample VCF containing transcript, Wildtype protein sequence, and Downstream protein sequence information",
        required=True,
    )
    fasta_parser.add_argument(
        "-l",
        "--length",
        action="store",
        dest="length",
        type=int,
        default=21,
        help="Length of the peptide sequence to use when creating the FASTA (default: %(default)s)",
    )
    fasta_parser.add_argument(
        "-d",
        "--downstream_length",
        action="store",
        dest="downstream_length",
        default=1000,
        help="Cap to limit the downstream sequence length for frameshifts when creating the fasta file."
        + "Use 'full' to include the full downstream sequence. (Default: %(default)s)",
    )
    fasta_parser.add_argument(
        "-o",
        "--output",
        action="store",
        dest="output",
        default="output.fasta",
        help="output fasta file name (default: %(default)s)",
    )

    return parser


def main():
    parser = parse_args()
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()

    external_tool_checking()
    if args.sub_command == "indel":
        bwa_bam = preprocessing(args.input, ref=args.ref)
        scansv_caller(bwa_bam, ref=args.ref)
    elif args.sub_command == "anno":
        out_vcf = args.output
        pre_vcf = f"pre_vep.{os.getpid()}.vcf"
        vcf_renewer(
            in_vcf=args.input, out_vcf=pre_vcf, ref=args.ref, slippage=args.slippage
        )
        vep_caller(in_vcf=pre_vcf, out_vcf=out_vcf, ref=args.ref, cutoff=args.cutoff)
        remove(pre_vcf)
    elif args.sub_command == "hla":
        flag = False
        if not args.alleles and not args.bam:
            sys.stderr.write("HLA class I alleles [OR] RNA-Seq BAM file is needed!")
            sys.exit(0)
        elif not args.alleles and args.bam:
            alleles = OptiType_runner(args.bam).split(",")
            flag = True
        elif args.alleles and not args.bam:
            alleles = args.alleles.split(",")
            flag = True
        else:
            sys.stderr.write(
                "HLA alleles [OR] RNA-seq BAM file is needed! Cannot process both infomation!"
            )
            sys.exit(0)
        if flag:
            epitope_lengths = list(map(int, args.epitope_lengths.split(",")))
            path_to_iedb = args.path_to_iedb
            peptide_sequence_length = args.length
            af_field = args.af_field
            vcf = args.vcf
            status_message(f"Begin to process VCF: {vcf}")
            status_message("......")
            start_execute_timestamp = datetime.datetime.now()
            if not args.output:
                out_filename = "{}.tsv".format(os.path.basename(vcf).split(".")[0])
            else:
                out_filename = args.output

            temp_dir = generate_protein_fasta(
                input_vcf=vcf,
                peptide_sequence_length=peptide_sequence_length,
                epitope_lengths=epitope_lengths,
                downstream_length=1000,
            )
            result_list = call_iedb_and_parse_outputs(
                path_to_iedb=path_to_iedb,
                alleles=alleles,
                epitope_lengths=epitope_lengths,
                temp_dir=temp_dir,
                methods=["ann", "netmhcpan"],
                top_score_metric=args.metric,
                number_of_threads=args.thread,
            )
            combined_filename = f"{os.getpid()}.combined.txt"
            combined_file = combine_results(
                infile_list=result_list,
                outfile=combined_filename,
                top_score_metric=args.metric,
                binding_cutoff=args.binding_threshold,
            )
            if combined_file:
                ScanNeo_utils.add_ranking(
                    infile=combined_file,
                    outfile=out_filename,
                    invcf=vcf,
                    af_field=af_field,
                )
                remove(combined_filename)
            end_execute_timestamp = datetime.datetime.now()
            elapsed_time = (
                end_execute_timestamp - start_execute_timestamp
            ).total_seconds()
            status_message(
                f"Processing {vcf} completed, consumed {elapsed_time} seconds"
            )
    elif args.sub_command == "fasta":
        generate_protein_fasta_main(
            input_vcf=args.vcf,
            peptide_sequence_length=args.length,
            output_fasta_file=args.output,
            downstream_sequence_length=args.downstream_length,
        )


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(1)
