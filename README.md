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
```$ sudo apt-get install csh gawk```

Install them on CentOS/RHEL.
```# yum install tcsh gawk```


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
	vcf_file			:VCF file with INDELs in it. (rna-seq.indel.vcf)
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
	
	input_bam_file   			:input BAM file is produced by transIndel_build.py
	
#### Output:
	
	output_vcf_file   			:Reported Indels with VCF format
  
#### STEP 3: neoantigen prediction
```
ScanNeo.py hla -i vep.vcf --alleles allele1,allele2 -e 8,9 -o output.tsv [options]
ScanNeo.py hla -b RNA_seq.bam --alleles allele1,allele2 -e 8,9 -o output.tsv [options]
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

#### Input:
	
	input_bam_file   			:input BAM file is produced by transIndel_build.py
	
#### Output:
	
	output_vcf_file   			:Reported Indels with VCF format  

