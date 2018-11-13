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
Change from
```python
this_dir = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(os.path.realpath(__file__))))
```
to
```python
this_dir = os.path.dirname(os.path.realpath(__file__))
```
##### configure config.ini of Optitype
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
#### configure yara
Index the HLA reference genome hla_reference_rna.fasta from Optitype
path/to/anaconda2/share/optitype-1.3.2-1/data/hla_reference_rna.fasta
```bash
yara_index hla_reference_rna.fasta -o hla.index  
```
By executing this command, it will generate 12 files, namely, 
hla.index.lf.drp
hla.index.lf.drv
hla.index.rid.concat
hla.index.sa.ind
hla.index.sa.val
hla.index.txt.limits
hla.index.lf.drs
hla.index.lf.pst
hla.index.rid.limits
hla.index.sa.len
hla.index.txt.concat
hla.index.txt.size
#### configure config.ini file
```
[fasta]

hg38=/path/to/hg38.fa
hg19=/path/to/hg19.fa

[annotation]

hg38=/path/to/gencode.v21.annotation.gtf
hg19=/path/to/gencode.v19.annotation.gtf

[yara]

index=/path/to/hla.index
```

