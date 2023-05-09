import pandas as pd
import os

# Fungi samples paired
fun_pd = pd.read_csv("raw/proj_no_literature_incomplete/fungi_paired.tsv", header = 0, names = ["run_id"], sep = "\t")
fun_pd_id = list(fun_pd["run_id"])
#fun_pd_id = fun_pd_id[0:2]

# Fungi samples single
fun_sg = pd.read_csv("raw/proj_no_literature_incomplete/fungi_single.tsv", header = 0, names = ["run_id"], sep = "\t")
fun_sg_id = list(fun_sg["run_id"])
#fun_sg_id = fun_sg_id[0:2]

# Cache directories
CACHE_FUN_FQ_PD = "01-fun-fq-pd"
CACHE_FUN_FQ_SG = "01-fun-fq-sg"

CACHE_FUN_MERGED = "02-fun-merged"
CACHE_PD_SG_2GETHER = "03-fastq-merged-sg"
CACHE_FQ2FA = "04-fq-2-fa"
CACHE_ITS = "05-extracted-its"

# Functions to generate file names
def get_fun_fq_pd_files(wildcards):
    global FQ
    FQ = expand(os.path.join("cache", CACHE_FUN_FQ_PD, "{ID}_{PD}.fastq"), ID = fun_pd_id, PD = [1,2])
    return FQ

def get_fun_fq_sg_files(wildcards):
    global FQ
    FQ = expand(os.path.join("cache", CACHE_FUN_FQ_SG, "{ID}.fastq"), ID = fun_sg_id)
    return FQ

def get_merged(wildcards):
    global FQ
    FQ = expand(os.path.join("cache", CACHE_FUN_MERGED, "{ID}.fastq"), ID = fun_pd_id)
    return FQ

def get_merged_sg_files(wildcards):
    ck_output = checkpoints.paste_pd_sg_together.get(**wildcards).output[0]
    global SMP
    SMP, = glob_wildcards(os.path.join(ck_output, "{SAMPLE}.fastq"))
    merged_sg_files = expand(os.path.join(ck_output, "{SAMPLE}.fastq"), SAMPLE=SMP)
    return merged_sg_files

def get_fasta(wildcards):
    ck_output = checkpoints.paste_pd_sg_together.get(**wildcards).output[0]
    global SMP
    SMP, = glob_wildcards(os.path.join(ck_output, "{SAMPLE}.fastq"))
    fa = expand(os.path.join("cache", CACHE_FQ2FA, "{SAMPLE}.fasta"), SAMPLE=SMP)
    return fa

def get_its_files(wildcards):
    ck_output = checkpoints.paste_pd_sg_together.get(**wildcards).output[0]
    global SMP
    SMP, = glob_wildcards(os.path.join(ck_output, "{SAMPLE}.fastq"))
    ITS = expand(os.path.join("cache", CACHE_ITS, "{SAMPLE}.{ITS}.fasta"), SAMPLE = SMP, ITS = ["ITS1", "ITS2"])
    return ITS

rule all:
    input:
        get_fun_fq_pd_files,
        get_fun_fq_sg_files,
        get_merged,
        get_fasta,
        get_its_files,
        "cache/finished.txt"
        
# Download paired end files
rule raw_fastq_pd:
    output:
        os.path.join("cache", CACHE_FUN_FQ_PD, "{ID}_1.fastq"),
        os.path.join("cache", CACHE_FUN_FQ_PD, "{ID}_2.fastq")
    conda:
        "src/conda/sra-tools.loose.yaml"
    threads: 1
    params:
        os.path.join("cache", CACHE_FUN_FQ_PD)
    shell:
        """
        fasterq-dump \
            --split-3 \
            --threads {threads} \
            -O {params} \
            -f \
            {wildcards.ID}
        """
# Download single end files
rule raw_fastq_sg:
    output:
        os.path.join("cache", CACHE_FUN_FQ_SG, "{ID}.fastq")
    conda:
        "src/conda/sra-tools.loose.yaml"
    threads: 1
    params:
        os.path.join("cache", CACHE_FUN_FQ_SG)
    shell:
        """
        fasterq-dump \
            --split-3 \
            --threads {threads} \
            -O {params} \
            {wildcards.ID}
        """

# Merge paired end reads
rule merge_pd:
    input:
        r1=os.path.join("cache", CACHE_FUN_FQ_PD, "{ID}_1.fastq"),
        r2=os.path.join("cache", CACHE_FUN_FQ_PD, "{ID}_2.fastq")
    output:
        os.path.join("cache", CACHE_FUN_MERGED, "{ID}.fastq")
    threads: 1
    singularity:
        "https://depot.galaxyproject.org/singularity/ngmerge:0.3--ha92aebf_1"
    shell:
        """
        NGmerge \
            -1 {input.r1} \
            -2 {input.r2} \
            -n {threads} \
            -o {output}
        """

# Paste merged paired and single end files together
checkpoint paste_pd_sg_together:
    input:
        get_merged,
        get_fun_fq_sg_files
    output:
        directory(os.path.join("cache", CACHE_PD_SG_2GETHER))
    threads: 1
    shell:
        """
        mkdir -p {output}
        cp {input} {output}
        """
        
# Convert fastq to fasta
rule fastq2fasta:
    input:
        os.path.join("cache", CACHE_PD_SG_2GETHER , "{SAMPLE}.fastq")
    output:
        os.path.join("cache", CACHE_FQ2FA, "{SAMPLE}.fasta")
    threads: 1
    shell:
        """
        sed -n '1~4s/^@/>/p;2~4p' {input} > {output}
        """

# Extract ITS regions
rule get_ITS1_ITS2:
    input:
        os.path.join("cache", CACHE_FQ2FA, "{SAMPLE}.fasta")
    output:
       os.path.join("cache", CACHE_ITS, "{SAMPLE}.ITS1.fasta"),
       os.path.join("cache", CACHE_ITS, "{SAMPLE}.ITS2.fasta")
    singularity:
        "docker://quay.io/biocontainers/itsx:1.1.3--hdfd78af_1"
    threads: 1
    resources:
        mem_mb=50000
    params:
        os.path.join("cache", CACHE_ITS)
    shell:
        """
        ITSx \
            --cpu {threads} \
            -i {input} \
            -o {params}/{wildcards.SAMPLE}
        """

#####################################################################################################################################
# 16S analysis

# Bacteria samples paired
bac_pd = pd.read_csv("raw/proj_no_literature_incomplete/bacteria_paired.tsv", header = 0, names = ["run_id"], sep = "\t")
bac_pd_id = list(bac_pd["run_id"])
#bac_pd_id = bac_pd_id[135:140]

# Bacteria samples single
bac_sg = pd.read_csv("raw/proj_no_literature_incomplete/bacteria_single.tsv", header = 0, names = ["run_id"], sep = "\t")
bac_sg_id = list(bac_sg["run_id"])
#bac_sg_id = bac_sg_id[55:60]

# Cache directories
CACHE_BAC_FQ_PD = "01-bac-fq-pd"
CACHE_BAC_FQ_SG = "01-bac-fq-sg"

CACHE_BAC_MERGED = "02-bac-merged"
CACHE_BAC_PD_SG_2GETHER = "03-bac-fastq-merged-sg"
CACHE_BAC_FQ2FA = "04-bac-fq-2-fa"
CACHE_16S = "05-aligned-16s"
CACHE_BWA_DB = "06-bwa-db"
CACHE_ALN = "07-bwa-aln"
CACHE_COV = "08-coverage"
CACHE_READ_LENGTH = "09-reads-length"
CACHE_BAM2BED = "10-bam2bed"


# Functions to generate file names
def get_bac_fq_pd_files(wildcards):
    global FQ
    FQ = expand(os.path.join("cache", CACHE_BAC_FQ_PD, "{ID}_{PD}.fastq"), ID = bac_pd_id, PD = [1,2])
    return FQ

def get_bac_fq_sg_files(wildcards):
    global FQ
    FQ = expand(os.path.join("cache", CACHE_BAC_FQ_SG, "{ID}.fastq"), ID = bac_sg_id)
    return FQ

def get_bac_merged(wildcards):
    global FQ
    FQ = expand(os.path.join("cache", CACHE_BAC_MERGED, "{ID}.fastq"), ID = bac_pd_id)
    return FQ

def get_bac_merged_sg_files(wildcards):
    ck_output = checkpoints.paste_bac_pd_sg_together.get(**wildcards).output[0]
    global SMP3
    SMP3, = glob_wildcards(os.path.join(ck_output, "{SAMPLE}.fastq"))
    merged_sg_files = expand(os.path.join(ck_output, "{SAMPLE}.fastq"), SAMPLE=SMP3)
    return merged_sg_files

def get_bwa_aln(wildcards):
    ck_output = checkpoints.paste_bac_pd_sg_together.get(**wildcards).output[0]
    global SMP
    SMP, = glob_wildcards(os.path.join(ck_output, "{SAMPLE}.fastq"))
    fa = expand(os.path.join("cache", CACHE_ALN, "{SAMPLE}.bam"), SAMPLE=SMP)
    return fa

def get_coverage(wildcards):
    ck_output = checkpoints.paste_bac_pd_sg_together.get(**wildcards).output[0]
    global SMP
    SMP, = glob_wildcards(os.path.join(ck_output, "{SAMPLE}.fastq"))
    fa = expand(os.path.join("cache", CACHE_COV, "{SAMPLE}.txt"), SAMPLE=SMP)
    return fa

def get_reads_length(wildcards):
    ck_output = checkpoints.paste_bac_pd_sg_together.get(**wildcards).output[0]
    global SMP
    SMP, = glob_wildcards(os.path.join(ck_output, "{SAMPLE}.fastq"))
    fq = expand(os.path.join("cache", CACHE_READ_LENGTH, "{SAMPLE}.txt"), SAMPLE=SMP)
    return fq

def get_bed(wildcards):
    ck_output = checkpoints.paste_bac_pd_sg_together.get(**wildcards).output[0]
    global SMP
    SMP, = glob_wildcards(os.path.join(ck_output, "{SAMPLE}.fastq"))
    bed = expand(os.path.join("cache", CACHE_BAM2BED, "{SAMPLE}.bed"), SAMPLE=SMP)
    return bed

# Download paired end files
rule raw_bac_fastq_pd:
    output:
        os.path.join("cache", CACHE_BAC_FQ_PD, "{ID}_1.fastq"),
        os.path.join("cache", CACHE_BAC_FQ_PD, "{ID}_2.fastq")
    conda:
        "src/conda/sra-tools.loose.yaml"
    threads: 1
    params:
        os.path.join("cache", CACHE_BAC_FQ_PD)
    shell:
        """
        fasterq-dump \
            --split-3 \
            --threads {threads} \
            -O {params} \
            -f \
            {wildcards.ID}
        """
#Download single end files
rule raw_bac_fastq_sg:
    output:
        os.path.join("cache", CACHE_BAC_FQ_SG, "{ID}.fastq")
    conda:
        "src/conda/sra-tools.loose.yaml"
    threads: 1
    params:
        os.path.join("cache", CACHE_BAC_FQ_SG)
    shell:
        """
        fasterq-dump \
            --split-3 \
            --threads {threads} \
            -O {params} \
            {wildcards.ID}
        """

# Merge paired end reads
rule merge_bac_pd:
    input:
        r1=os.path.join("cache", CACHE_BAC_FQ_PD, "{ID}_1.fastq"),
        r2=os.path.join("cache", CACHE_BAC_FQ_PD, "{ID}_2.fastq")
    output:
        os.path.join("cache", CACHE_BAC_MERGED, "{ID}.fastq")
    threads: 1
    params:
        "raw/qual_score.txt"
    singularity:
        "docker://staphb/bbtools:39.01"
    shell:
        """
        bbmerge.sh \
            in1={input.r1} \
            in2={input.r2} \
            out={output}
        """

# Paste merged paired and single end files together
checkpoint paste_bac_pd_sg_together:
    input:
        get_bac_merged,
        get_bac_fq_sg_files
    output:
        directory(os.path.join("cache", CACHE_BAC_PD_SG_2GETHER))
    threads: 1
    shell:
        """
        mkdir -p {output}
        cp {input} {output}
        """

rule create_aligner_db:
    input:
        "ref/silva138.1ssuRef_bacteria_16sHypervariableRegions_headerfix.fa"
    output:
        multiext(os.path.join("cache", CACHE_BWA_DB, "16hyperVar"),
            ".pac",
            ".ann",
            ".amb",
            ".bwt",
            ".sa")
    threads: 1
    log:
        os.path.join("cache/logs", CACHE_BWA_DB, "createDB.log")
    params:
        os.path.join("cache", CACHE_BWA_DB, "16hyperVar")
    singularity:
        "https://depot.galaxyproject.org/singularity/bwakit:0.7.17.dev1--hdfd78af_1"
    shell:
        """
        bwa index -p {params} {input}
        """

rule detect_16s_primer:
    input:
        fq=os.path.join("cache", CACHE_BAC_PD_SG_2GETHER , "{SAMPLE}.fastq"),
        db=rules.create_aligner_db.output
    output:
        os.path.join("cache", CACHE_ALN, "{SAMPLE}.bam")
    threads: 12
    resources:
        mem_mb=50000
    params:
        db=rules.create_aligner_db.params
    singularity:
        "https://depot.galaxyproject.org/singularity/bwakit:0.7.17.dev1--hdfd78af_1"
    shell:
        """
        bwa mem -t {threads} {params.db} {input.fq} | samtools sort -@{threads} -o {output} -
        """

rule coverage_16s:
    input:
        os.path.join("cache", CACHE_ALN, "{SAMPLE}.bam")
    output:
        os.path.join("cache", CACHE_COV, "{SAMPLE}.txt")
    threads: 24
    resources:
        mem_mb=50000
    singularity:
        "docker://staphb/bbtools:39.01"
    shell:
        """
        pileup.sh \
            nzo=t \
            in={input} \
            out={output}
        """

rule read_length:
    input:
        os.path.join("cache", CACHE_BAC_PD_SG_2GETHER , "{SAMPLE}.fastq")
    output:
        os.path.join("cache", CACHE_READ_LENGTH, "{SAMPLE}.txt")
    threads: 1
    singularity:
        "docker://quay.io/biocontainers/seqkit:2.3.1--h9ee0642_0"
    shell:
        """
        seqkit fx2tab -nl {input} > {output}
        """

rule bam2bed:
    input:
        os.path.join("cache", CACHE_ALN, "{SAMPLE}.bam")
    output:
        os.path.join("cache", CACHE_BAM2BED, "{SAMPLE}.bed")
    threads: 12
    resources:
        mem_mb=50000
    singularity:
        "https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h468198e_3"
    shell:
        """
        bedtools bamtobed -i {input} > {output}
        """

rule calculate_read_coverage:
    input:
        get_coverage,
        get_reads_length,
        get_bed
    output:
        "logs/calculate_read_coverage.log"
    threads: 12
    resources:
        mem_mb=100000
    singularity:
        "docker://rocker/tidyverse:4.2.2"
    params:
        "src/0r-microbiome/analyze.R"
    shell:
        """
        Rscript {params} > {output}
        """

rule done_16s:
    input:
        get_bac_fq_pd_files,
        get_bac_fq_sg_files,
        get_bac_merged,
        #get_bac_fasta,
        get_bac_merged_sg_files,
        get_bwa_aln,
        get_coverage,
        get_reads_length,
        get_bed,
        "logs/calculate_read_coverage.log"
    output:
        "cache/finished.txt"
    shell:
        """
        touch {output}
        """
