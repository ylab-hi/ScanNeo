Introduction
------------
a pipeline for personalized neoantigen discovery using RNA-seq data

Getting Started
----------------
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

Prerequisites
----------------
You need Python 2.7 to run ScanNeo.

### install necessary python packages via anaconda
Install [anaconda](https://www.anaconda.com/download/) firstly, then install dependent packages via conda in bioconda channel.
```
conda install -c bioconda optitype
conda install -c bioconda ensembl-vep
conda install -c bioconda sambamba
conda install -c bioconda bedtools
conda install -c bioconda picard
conda install -c bioconda bwa
conda install -c bioconda yara
conda install -c bioconda razers3
conda install -c conda-forge glpk
```
### install transIndel
```
git clone https://github.com/cauyrd/transIndel 
```

Place __transIndel_build_RNA.py__ and __transIndel_call.py__ in your folder in PATH.

### intall IEDB binding prediction tools
Download the archives for [HLA class I](http://tools.iedb.org/mhci/download/) and unpack them
```
wget -c https://downloads.iedb.org/tools/mhci/2.19.1/IEDB_MHC_I-2.19.1.tar.gz
tar -zxvf IEDB_MHC_I-2.19.1.tar.gz
cd mhc_i
./configure
```

__tcsh__ and __gawk__ are needed for IEDB binding prediction tools. You have to install them if they are not available.

Install them on Debian/Ubuntu/Mint Linux.
```$ sudo apt-get install tcsh gawk```

Install them on CentOS/RHEL.
```# yum install tcsh gawk```

### install VEP plugins

```
git clone https://github.com/ylab-hi/ScanNeo.git
cd VEP_plugins
cp Downstream.pm ~/.vep/Plugins
cp Wildtype.pm ~/.vep/Plugins
```


Configuration
----------------
### configure Optitype

Make modification on __OptiTypePipeline.py__
Change from
```python
this_dir = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(os.path.realpath(__file__))))
```
to
```python
this_dir = os.path.dirname(os.path.realpath(__file__))
```
#### configure config.ini of Optitype
```
[mapping]

# Absolute path to RazerS3 binary, and number of threads to use for mapping
razers3=/path/to/anaconda2/bin/razers3
threads=16

[ilp]
# A Pyomo-supported ILP solver. The solver must be globally accessible in the
# environment OptiType is run, so make sure to include it in PATH.
# Note: this is NOT a path to the solver binary, but a keyword argument for
# Pyomo. Examples: glpk, cplex, cbc.

solver=glpk
threads=1
```
### configure yara
Index the HLA reference genome hla_reference_rna.fasta from Optitype
path/to/anaconda2/share/optitype-1.3.2-1/data/hla_reference_rna.fasta
```bash
yara_index hla_reference_rna.fasta -o hla.index  
```
By executing this command, it will generate 12 files, namely, 
* hla.index.lf.drp
* hla.index.lf.drv
* hla.index.rid.concat
* hla.index.sa.ind
* hla.index.sa.val
* hla.index.txt.limits
* hla.index.lf.drs
* hla.index.lf.pst
* hla.index.rid.limits
* hla.index.sa.len
* hla.index.txt.concat
* hla.index.txt.size
### configure config.ini file
```
[fasta]

# reference genome file in FASTA format

hg38=/path/to/hg38.fa
hg19=/path/to/hg19.fa

[annotation]

# gene annotation file in GTF format

hg38=/path/to/gencode.v21.annotation.gtf
hg19=/path/to/gencode.v19.annotation.gtf

[yara]

# yara index for 

index=/path/to/hla.index
```

Usage
-------------------------
#### STEP 1: INDEL calling using RNA-seq data
```
ScanNeo.py indel -i rnaseq_bam -r hg38
```

#### Options:

```	
-h, --help            show this help message and exit
-i INPUT, --input INPUT
                        RNA-seq alignment file (BAM)
-r {hg19,hg38}, --ref {hg19,hg38}
                        reference genome (default: hg38)
```

#### Input:
```	
input_bam_file      :input RNA-seq BAM file. (e.g., rna-seq.bam)
reference_genome    :specify reference genome (hg19 or hg38)
```

#### Output:
```
	
your_output_bam_file		:BAM file for CIGAR string redefinement. (rna-seq.indel.bam)
vcf_file			:Reported Indels with VCF format. (rna-seq.indel.vcf)
```

#### STEP 2: indels annotation using VEP
```
ScanNeo.py anno -i input_vcf_file -o output_annotated_vcf_file [options]	
```

#### Options:
```
-h, --help            show this help message and exit
-i INPUT, --input INPUT
                        input VCF file
-c CUTOFF, --cutoff CUTOFF
                        MAF cutoff default: 0.01
-r {hg19,hg38}, --ref {hg19,hg38}
                        reference genome (default: hg38)
-o OUTPUT, --output OUTPUT
                        output annotated and filtered vcf file (default:
                        output.vcf)
```

#### Input:
```	
input_vcf_file   	        :input VCF file is produced by ScanNeo indel
cutoff   			:MAF cutoff according to 1000 genome project and gnomAD project
reference_genome                :specify reference genome (hg19 or hg38)
```

#### Output:

```	
output_annotated_vcf_file   			:VEP annotated Indels with VCF format
```

#### STEP 3: neoantigen prediction
```
ScanNeo.py hla -i vep.vcf --alleles HLA-A*02:01,HLA-B*08:01,HLA-C*03:03 -e 8,9 -o output.tsv [options]
ScanNeo.py hla -i vep.vcf -b RNA_seq.bam -e 8,9 -o output.tsv [options]
```

#### Options:
```
-h, --help            show this help message and exit
-i VCF, --input VCF   VEP annotated and filtered VCF
--alleles ALLELES     Name of the allele to use for epitope prediction.
                        Multiple alleles can be specified using a comma-
                        separated listinput HLA class I alleles
-b BAM, --bam BAM     Input RNA-Seq BAM file if you don't know sample HLA class I alleles
-l LENGTH, --length LENGTH
                         Length of the peptide sequence to use when creating
                         the FASTA (default: 21)
--binding BINDING_THRESHOLD
                         binding threshold ic50 (default: 500 nM)
-e EPITOPE_LENGTHS, --epitope-length EPITOPE_LENGTHS
                         Length of subpeptides (neoepitopes) to predict.
                         Multiple epitope lengths can be specified using a
                         comma-separated list. Typical epitope lengths vary
                         between 8-11. (default: 8,9,10,11)
-p PATH_TO_IEDB, --path-to-iedb PATH_TO_IEDB
                         Directory that contains the local installation of IEDB
                             
-m {lowest,median}, --metric {lowest,median}
                        The ic50 scoring metric to use when filtering epitopes
                        by binding-threshold lowest: Best MT Score - lowest MT
                        ic50 binding score of all chosen prediction methods.
                        median: Median MT Score - median MT ic50 binding score
                        of all chosen prediction methods. (default: lowest)
-o OUTPUT, --output OUTPUT
                       output text file name, name.tsv
```

#### Output:

```
name.tsv file contains neoantigen prediction results
```

__Report Columns__

|Column Name | Description |
| ---------- | ----------- |
|Chromosome  | The chromosome of this variant|
|Start       | The start position of this variant in the zero-based, half-open coordinate system |
|Stop	     | The stop position of this variant in the zero-based, half-open coordinate system |
|Reference   | The reference allele |
|Variant     | The alternate allele |
|Transcript  | The Ensembl ID of the affected transcript |
|Ensembl Gene ID |	The Ensembl ID of the affected gene |
|Variant Type |	The type of variant. missense for missense mutations, inframe_ins for inframe insertions, inframe_del for inframe deletions, and FS for frameshift variants |
|Mutation     |	The amnio acid change of this mutation |
|Protein Position |	The protein position of the mutation |
|Gene Name    | The Ensembl gene name of the affected gene |
|HLA Allele   | The HLA allele for this prediction |
|Peptide Length | The peptide length of the epitope |
|Sub-peptide Position | The one-based position of the epitope in the protein sequence used to make the prediction |
|Mutation Position    | The one-based position of the start of the mutation in the epitope. 0 if the start of the mutation is before the epitope |
|MT Epitope Seq       | Mutant epitope sequence |
|WT Epitope Seq   | Wildtype (reference) epitope sequence at the same position in the full protein sequence. NA if there is no wildtype sequence at this position or if more than half of the amino acids of the mutant epitope are mutated |
|Best MT Score Method | Prediction algorithm with the lowest mutant ic50 binding affinity for this epitope |
|Best MT Score        | Lowest ic50 binding affinity of all prediction algorithms used | 
|Corresponding WT Score | ic50 binding affinity of the wildtype epitope. NA if there is no WT Epitope Seq.|
|Corresponding Fold Change | Corresponding WT Score / Best MT Score. NA if there is no WT Epitope Seq.|
|Median MT Score | Median ic50 binding affinity of the mutant epitope of all prediction algorithms used |
|Median WT Score | Median ic50 binding affinity of the wildtype epitope of all prediction algorithms used. NA if there is no WT Epitope Seq. |
|Median Fold Change | Median WT Score / Median MT Score. NA if there is no WT Epitope Seq. |
|Individual Prediction Algorithm WT and MT Scores (multiple) | ic50 scores for the MT Epitope Seq and WT Eptiope Seq for the individual prediction algorithms used |

