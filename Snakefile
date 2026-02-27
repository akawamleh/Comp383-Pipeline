#pipelineproject
#snakemake workflow for steps2-5 with test/full mode
#final output is PipelineReport.txt
#auto-downloads reference and builds blast db if needed

import os

configfile: "my_config.yaml"

MODE=config["mode"]
SAMPLES=config["samples"]
THREADS=int(config["threads"])
KMER=int(config["spades_k"])

TESTDIR=config["paths"]["test_dir"]
FULLDIR=config["paths"]["full_dir"]

HCMV_ACC=config["reference"]["hcmv_accession"]
HCMV_FNA=config["reference"]["hcmv_fna"]
HCMV_INDEX=config["reference"]["bowtie2_index_prefix"]

BETA_TAXON=config["blast"]["taxon"]
BETA_FNA=config["blast"]["fasta"]
BLAST_PREFIX=config["blast"]["db_prefix"]

RESULTS=f"results/{MODE}"

#choose input directory based on mode
def raw_dir():
    return TESTDIR if MODE=="test" else FULLDIR

def fastq_path(sample,mate):
    d=raw_dir()
    return f"{d}/{sample}_{mate}.fastq.gz"

rule all:
    input:
        "PipelineReport.txt"

#refdownload

#download hcmv reference genome fasta using ncbi datasets
rule fetch_hcmv_ref:
    output:
        HCMV_FNA
    params:
        acc=HCMV_ACC
    shell:
        r"""
        mkdir -p ref
        rm -rf ref/hcmv_dataset
        rm -f ref/hcmv.zip

        datasets download genome accession {params.acc} --include genome -o ref/hcmv.zip
        unzip -o ref/hcmv.zip -d ref/hcmv_dataset >/dev/null

        fna=$(find ref/hcmv_dataset -name "*genomic.fna" | head -n 1)
        if [ -z "$fna" ]; then
          echo "ERROR:no genomic.fna found for HCMV" >&2
          exit 1
        fi

        cp "$fna" {output}
        """

#build bowtie2 index for hcmv reference
rule bowtie2_index:
    input:
        HCMV_FNA
    output:
        HCMV_INDEX + ".1.bt2",
        HCMV_INDEX + ".2.bt2",
        HCMV_INDEX + ".3.bt2",
        HCMV_INDEX + ".4.bt2",
        HCMV_INDEX + ".rev.1.bt2",
        HCMV_INDEX + ".rev.2.bt2"
    shell:
        r"""
        bowtie2-build {input} {HCMV_INDEX}
        """

#blastdbbuild

#download betaherpesvirinae genomes fasta using ncbi datasets
rule fetch_betaherpes_fna:
    output:
        BETA_FNA
    params:
        taxon=BETA_TAXON
    shell:
        r"""
        mkdir -p blast_db
        rm -rf blast_db/beta_dataset
        rm -f blast_db/beta.zip

        datasets download virus genome taxon {params.taxon} --refseq --include genome -o blast_db/beta.zip
        unzip -o blast_db/beta.zip -d blast_db/beta_dataset >/dev/null

        #datasets output fasta path
        src="blast_db/beta_dataset/ncbi_dataset/data/genomic.fna"
        if [ ! -s "$src" ]; then
          echo "ERROR:could not find or empty genomic.fna at $src" >&2
          exit 1
        fi

        cp "$src" {output}
        """

#make nucleotide blast database
rule make_blast_db:
    input:
        BETA_FNA
    output:
        BLAST_PREFIX + ".nsq"
    shell:
        r"""
        makeblastdb -in {input} -out {BLAST_PREFIX} -title betaherpesvirinae -dbtype nucl
        """

#step2mapandfilter

#map reads to hcmv reference and write sorted bam
rule map_sorted_bam:
    input:
        idx=rules.bowtie2_index.output,
        r1=lambda wc: fastq_path(wc.sample,1),
        r2=lambda wc: fastq_path(wc.sample,2)
    output:
        bam=f"{RESULTS}/bam/{{sample}}.sorted.bam"
    threads: THREADS
    shell:
        r"""
        mkdir -p {RESULTS}/bam
        bowtie2 -p {threads} -x {HCMV_INDEX} -1 {input.r1} -2 {input.r2} |
        samtools view -bS - |
        samtools sort -@ {threads} -o {output.bam}
        """

#filter to properly paired mapped reads and export paired fastq.gz
rule filter_to_mapped_fastq:
    input:
        bam=f"{RESULTS}/bam/{{sample}}.sorted.bam"
    output:
        r1=f"{RESULTS}/mapped_fastq/{{sample}}_1.fastq.gz",
        r2=f"{RESULTS}/mapped_fastq/{{sample}}_2.fastq.gz"
    threads: THREADS
    shell:
        r"""
        mkdir -p {RESULTS}/mapped_fastq

        samtools view -b -f 2 -F 4 {input.bam} > {RESULTS}/bam/{wildcards.sample}.mapped_proper.bam
        samtools sort -n -@ {threads} -o {RESULTS}/bam/{wildcards.sample}.mapped_proper.namesorted.bam {RESULTS}/bam/{wildcards.sample}.mapped_proper.bam

        samtools fastq -@ {threads} \
          -1 {output.r1} -2 {output.r2} \
          -0 /dev/null -s /dev/null -n \
          {RESULTS}/bam/{wildcards.sample}.mapped_proper.namesorted.bam
        """

#step2counts

#count pairs before filtering (from R1 input)
rule count_before:
    input:
        r1=lambda wc: fastq_path(wc.sample,1)
    output:
        f"{RESULTS}/counts/{{sample}}_before.txt"
    shell:
        r"""
        mkdir -p {RESULTS}/counts
        python scripts/count_pairs.py -i {input.r1} -o {output}
        """

#count pairs after filtering (from mapped R1)
rule count_after:
    input:
        r1=f"{RESULTS}/mapped_fastq/{{sample}}_1.fastq.gz"
    output:
        f"{RESULTS}/counts/{{sample}}_after.txt"
    shell:
        r"""
        mkdir -p {RESULTS}/counts
        python scripts/count_pairs.py -i {input.r1} -o {output}
        """

#step3assembly

#assemble filtered reads using spades
rule spades:
    input:
        r1=f"{RESULTS}/mapped_fastq/{{sample}}_1.fastq.gz",
        r2=f"{RESULTS}/mapped_fastq/{{sample}}_2.fastq.gz"
    output:
        contigs=f"{RESULTS}/assemblies/{{sample}}/contigs.fasta"
    threads: THREADS
    shell:
        r"""
        mkdir -p {RESULTS}/assemblies/{wildcards.sample}
        spades.py -1 {input.r1} -2 {input.r2} -k {KMER} -t {threads} -o {RESULTS}/assemblies/{wildcards.sample}
        test -s {output.contigs}
        """

#step4stats

#compute stats for contigs >1000bp
rule contig_stats:
    input:
        contigs=f"{RESULTS}/assemblies/{{sample}}/contigs.fasta"
    output:
        stats=f"{RESULTS}/assemblies/{{sample}}/contigs_over_1000_stats.txt"
    shell:
        r"""
        python scripts/contig_stats_over_1000.py -i {input.contigs} -o {output.stats}
        """

#step5longestandblast

#extract longest contig
rule longest_contig:
    input:
        contigs=f"{RESULTS}/assemblies/{{sample}}/contigs.fasta"
    output:
        f"{RESULTS}/assemblies/{{sample}}/longest_contig.fasta"
    shell:
        r"""
        python scripts/longest_contig.py -i {input.contigs} -o {output}
        """

#blast longest contig against betaherpesvirinae db
rule blastn_longest:
    input:
        query=f"{RESULTS}/assemblies/{{sample}}/longest_contig.fasta",
        db=rules.make_blast_db.output
    output:
        f"{RESULTS}/assemblies/{{sample}}/{{sample}}_blast.tsv"
    shell:
        r"""
        blastn \
          -query {input.query} \
          -db {BLAST_PREFIX} \
          -out {output} \
          -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle" \
          -max_hsps 1
        """

#top5 blast hits
rule blast_top5:
    input:
        f"{RESULTS}/assemblies/{{sample}}/{{sample}}_blast.tsv"
    output:
        f"{RESULTS}/assemblies/{{sample}}/{{sample}}_blast_top5.tsv"
    shell:
        r"""
        python scripts/blast_top5.py -i {input} -o {output}
        """

#finalreport

#write PipelineReport.txt
rule report:
    input:
        expand(f"{RESULTS}/counts/{{sample}}_before.txt",sample=SAMPLES),
        expand(f"{RESULTS}/counts/{{sample}}_after.txt",sample=SAMPLES),
        expand(f"{RESULTS}/assemblies/{{sample}}/contigs_over_1000_stats.txt",sample=SAMPLES),
        expand(f"{RESULTS}/assemblies/{{sample}}/{{sample}}_blast_top5.tsv",sample=SAMPLES)
    output:
        "PipelineReport.txt"
    params:
        samples=" ".join(SAMPLES),
        results=RESULTS
    shell:
        r"""
        python scripts/write_report_final.py -o {output} -s {params.samples} -r {params.results}
        """