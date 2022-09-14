# APPRIS RNA-Seq pipeline
- [APPRIS RNA-Seq pipeline](#appris-rna-seq-pipeline)
  - [Introduction](#introduction)
  - [Installation](#installation)
  - [Usage](#usage)
  - [Directory structure](#directory-structure)
  - [Author information and license](#author-information-and-license)
  - [Release History](#release-history)
  - [Contributing](#contributing)

## Introduction

The [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline presented here allow to align a batch of RNA-seq paired-end samples accounting for genomic features.

It has been thought to be runned in a High-performance computing environment, but it could be adapted depending on the computer capability.

It allows both [GENCODE](ftp://ftp.ebi.ac.uk/pub/databases/gencode) and [RefSeq](https://ftp.ncbi.nlm.nih.gov/refseq/) genome reference annotations.

Sequentially it executes:

- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/): To find and removes adapters. It uses the `cutadapt_pe` mode.
- [STAR](https://github.com/alexdobin/STAR): Spliced and referenced transcripts aligner.
- [Samtools](http://www.htslib.org/): To interact with the sequencing data.
- [featureCounts](http://subread.sourceforge.net/): To count reads to genomic features such as genes or exons.

*[QSplice](https://gitlab.com/fpozoc/qsplice) could also be runned after the pipeline in order to quantify the splice junctions coverage per transcript.*

## Installation

Run the silent installation of Miniconda in case you don't have this software in your Linux Environment

```sh
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
```

Once you have installed Miniconda/Anaconda, create a Python 3.7 environment. Then, install snakemake in your conda environment:

```sh
conda install -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

If [Slurm](https://slurm.schedmd.com/documentation.html) has been installed in the computing environment, install with:

```sh
sbatch -o log.txt -e err.txt -J smk-install -c 2 --mem=2G -t 200 --wrap "mamba create -c conda-forge -c bioconda -y -n snakemake snakemake"
```

## Usage

To execute the pipeline, first the user must prepare their samples.

```sh
git clone git@gitlab.com:fpozoc/appris_rnaseq.git
cd appris_rnaseq
```

Then, edit config and workflow as needed:

```sh
vim config/config.yaml
```

sample `config.yaml`, in which user must decide:

- To use [GENCODE](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/) or [RefSeq](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz) annotation GTF files.
- Which is the appropriate reference genome for this analysis. Please, read the [Heng Li reference](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) and this [Biostar answer](https://www.biostars.org/p/342482/) before taking a final decision.
- Which `cutadapt` and `STAR` parameters desire to select
- Which RNA-seq samples user wants to align. In this case we are using [E-MTAB-2836](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2836/).

```sh
annotations:
    GRCh38:
        g34:
            url: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.primary_assembly.annotation.gtf.gz # GENCODE 34 GRCh38 gtf annotation file
            enabled: True # Enabled to be runned
        rs109:
            url: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf.gz  
            enabled: True
    GRCh37:
        g19:
            url: XXXX # GENCODE 19 GRCh37 gtf annotation file
            enabled: False
genomes:
    GRCh38:
        url: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz 
    GRCh37: 
        url: ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
params:
    cutadapt_pe: "option1"
    STAR: "option2"
samples:
     E-MTAB-XXXX: # proyect
         ERRXXXXX: # sample
             "1": # replicate
                 r1: ERR315325_1.fastq.gz # FastQ read 1
                 r2: ERR315325_2.fastq.gz # FastQ read 2
```

Finally, execute workflow, deploy software dependencies via conda:

```sh
snakemake -n --use-conda
```

If [Slurm](https://slurm.schedmd.com/documentation.html) has been installed in the computing environment:

```sh
snakemake -n --use-conda --profile slurm
```

## Directory structure

```sh
├── config
│   └── config.yaml
├── envs
│   └── environment.yaml
├── log
├── out
├── README.md
├── seq
└── Snakefile
```

## Author information and license

Fernando Pozo ([@fpozoca](https://twitter.com/fpozoca) – fpozoc@cnio.es) and Tomás Di Domenico.

Project initially forked from [here](https://gitlab.com/bu_cnio/appris_rnaseq).

## Release History

* 1.0.0.

## Contributing

1. Fork it (<https://gitlab.com/fpozoc/appris_rnaseq.git>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request
