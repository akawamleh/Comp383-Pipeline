# COMP 383 Pipeline Project

Adam Kawamleh

## Project Overview

This project implements a fully reproducible bioinformatics workflow using Snakemake to analyze Human Cytomegalovirus (HCMV) sequencing data.

For each sample, the pipeline performs:

- Download of the HCMV reference genome from NCBI
- Construction of a Bowtie2 index
- Alignment of paired-end reads to the HCMV reference
- Filtering to retain only properly mapped read pairs
- Counting read pairs before and after filtering
- Genome assembly using SPAdes (k=99)
- Contig statistics calculation for contigs >1000 bp
- Extraction of the longest contig
- BLAST of the longest contig against a Betaherpesvirinae database
- Automated report generation

All steps are automated and reproducible through Snakemake.

## Data Source (NCBI)

The sequencing data used in this project was obtained from the NCBI Sequence Read Archive (SRA).

The following accession numbers were used:

- SRR5660030  
- SRR5660033  
- SRR5660044  
- SRR5660045  

Reference genome:

- Human Cytomegalovirus (HCMV)  
- Accession: GCF_000845245.1  

BLAST database:

- Betaherpesvirinae genomes downloaded using the NCBI datasets CLI 

## Dependencies

The following software must be installed and available in your PATH:

- snakemake
- bowtie2
- samtools
- spades
- blastn (BLAST+)
- makeblastdb
- NCBI SRA Toolkit (fasterq-dump)
- NCBI datasets CLI
- gzip
- unzip
- python

## Installation

Create and activate a conda environment:

conda create -n cmv_pipeline snakemake bowtie2 samtools spades blast ncbi-datasets-cli sra-tools unzip -c bioconda -c conda-forge  
conda activate cmv_pipeline  

## Downloading FASTQ Data (REQUIRED)

The raw sequencing data comes from NCBI SRA.

Create the raw data directory:

mkdir -p data/raw

Download the FASTQ files using `fasterq-dump`:

fasterq-dump SRR5660030 --split-files  
fasterq-dump SRR5660033 --split-files  
fasterq-dump SRR5660044 --split-files  
fasterq-dump SRR5660045 --split-files  

Move the files into the raw directory:

mv SRR5660030_*.fastq data/raw/  
mv SRR5660033_*.fastq data/raw/  
mv SRR5660044_*.fastq data/raw/  
mv SRR5660045_*.fastq data/raw/  

Compress the FASTQ files:

gzip data/raw/*.fastq  

Final files should look like:

data/raw/SRR5660030_1.fastq.gz  
data/raw/SRR5660030_2.fastq.gz  
data/raw/SRR5660033_1.fastq.gz  
data/raw/SRR5660033_2.fastq.gz  
data/raw/SRR5660044_1.fastq.gz  
data/raw/SRR5660044_2.fastq.gz  
data/raw/SRR5660045_1.fastq.gz  
data/raw/SRR5660045_2.fastq.gz  

## Test Data (for fast grading)

Create smaller FASTQ files:

mkdir -p data/test

zcat data/raw/SRR5660030_1.fastq.gz | head -n 40000 | gzip > data/test/SRR5660030_1.fastq.gz  
zcat data/raw/SRR5660030_2.fastq.gz | head -n 40000 | gzip > data/test/SRR5660030_2.fastq.gz  

Repeat for all samples.



### Test Mode (recommended for grading)

Edit config.yaml:

mode: test  

Run:

snakemake --cores 4  

### Full Mode (all data)

Edit config.yaml:

mode: full  

Run:

snakemake --cores 4  

## Output Description

The final output file:

PipelineReport.txt

Each report includes:

- Read pairs before and after Bowtie2 filtering
- Number of contigs >1000 bp
- Total base pairs in those contigs
- Top 5 BLAST hits

BLAST output columns:

sacc pident length qstart qend sstart send bitscore evalue stitle

## Reproducibility

After cloning this repository and installing dependencies, the entire workflow can be executed with:

snakemake --cores 4  

The pipeline automatically:
- downloads the HCMV reference genome
- builds the Bowtie2 index
- builds the BLAST database

## Notes

- SPAdes is run with k=99 as required
- Only properly paired mapped reads are retained
- Test mode is recommended for quick validation