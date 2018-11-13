Introduction
------------
a pipeline for personalized neoantigen discovery using RNA-seq data

Getting Started
----------------
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

Prerequisites
----------------

What things you need to install the software and how to install them

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

Configuration
----------------
#### configure Optitype

Make modification on OptiTypePipeline.py 
```python
this_dir = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(os.path.realpath(__file__))))
```
```python
this_dir = os.path.dirname(os.path.realpath(__file__))
```

#### configure yara
path_to_anaconda2/share/optitype-1.3.2-1/data/hla_reference_rna.fasta
```bash
yara_index hla_reference_rna.fasta -o hla.index  
```
#### configure config.ini file
```
[fasta]

hg38=/home/tywang/database/genome/hg38.fa
hg19=/home/tywang/database/genome/hg19.fa

[annotation]

hg38=/home/tywang/database/gencode/gencode.v21.annotation.gtf
hg19=/home/tywang/database/gencode/gencode.v19.annotation.gtf

[yara]

index=/home/tywang/database/genome/hla.index
```

