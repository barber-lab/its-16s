import pandas as pd
import itertools

#### FUNCTIONS

# Functions to generate file names
def get_genome_meta_files(wildcards):
    global FILES
    FILES = expand(os.path.join("cache", CACHE_META, "{GENUS}.jsonl"), GENUS = genus)
    return FILES

# Functions to generate file names
def get_genome_meta_files_specialists(wildcards):
    global FILES
    FILES = expand(os.path.join("cache", CACHE_SPECIAL, "{SPECIAL}.jsonl"), SPECIAL = special)
    return FILES

# Functions to generate file names
def get_genome_from_ncbi(wildcards):
    global FILES1, FILES2
    FILES1 = expand(os.path.join("cache", CACHE_GENOME, "{ID}_{KING}_{PRED}.fna"), ID = b_genome, KING = ["bacteria"], PRED = ["prodigal"])
    FILES2 = expand(os.path.join("cache", CACHE_GENOME, "{ID}_{KING}_{PRED}.fna"), ID = f_genome, KING = ["fungi"], PRED = ["glimmerhmm"])
    FILES3 = expand(os.path.join("cache", CACHE_GFF3, "{ID}_{KING}_{PRED}.gff"), ID = b_genome, KING = ["bacteria"], PRED = ["prodigal"])
    FILES4 = expand(os.path.join("cache", CACHE_GFF3, "{ID}_{KING}_{PRED}.gff"), ID = f_genome, KING = ["fungi"], PRED = ["glimmerhmm"])
    return FILES1 + FILES2 + FILES3 + FILES4

    # FILES1 = expand(os.path.join("cache", CACHE_GENOME, "{ID}_{KING}.fna"), ID = b_genome, KING = ["Bacteria"])

# Functions to generate file names
def get_BGCs(wildcards):
    global BGC1, BGC2
    BGC1 = expand(os.path.join("cache", CACHE_BGC, "{ID}/{ID}_{KING}_{PRED}.json"), ID = b_genome, KING = ["bacteria"], PRED = ["prodigal"])
    BGC2 = expand(os.path.join("cache", CACHE_BGC, "{ID}/{ID}_{KING}_{PRED}.json"), ID = f_genome, KING = ["fungi"], PRED = ["glimmerhmm"])
    return BGC1 + BGC2

# Functions to generate file names
def get_BGCs_gff(wildcards):
    global BGC1, BGC2
    BGC1 = expand(os.path.join("cache", CACHE_BGC2, "{ID}/{ID}_{KING}_{PRED}.json"), ID = b_genome, KING = ["bacteria"], PRED = ["prodigal"])
    BGC2 = expand(os.path.join("cache", CACHE_BGC2, "{ID}/{ID}_{KING}_{PRED}.json"), ID = f_genome, KING = ["fungi"], PRED = ["glimmerhmm"])
    return BGC1 + BGC2

# Functions to generate file names
def get_protein(wildcards):
    global PROT
    PROT = expand(os.path.join("cache", CACHE_PROTEIN, "{ID}_{KING}_{PRED}.fna"), ID = b_genome, KING = ["bacteria"], PRED = ["prodigal"])
    return PROT

# Functions to generate file names
def get_amr(wildcards):
    global AMR
    AMR = expand(os.path.join("cache", CACHE_AMR, "{ID}_bacteria_{PRED}_amr.tsv"), ID = b_genome, PRED = ["prodigal"])
    return AMR


# def get_genome_names(wildcards):
#     ck_output = rules.download_genome.params[0]
#     global SMP
#     SMP, = glob_wildcards(os.path.join(ck_output, "{GENOME}.fna"))
#     FILES = expand(os.path.join(ck_output, "{SAMPLE}.fna"), SAMPLE=SMP)
#     return FILES

# def get_BGCs(wildcards):
#     ck_output = rules.download_genome.params[0]
#     global SMP
#     SMP, = glob_wildcards(os.path.join(ck_output, "{GENOME}.fna"))
#     FILES = expand(os.path.join("cache", CACHE_BGC, "{SAMPLE}.res"), SAMPLE=SMP)
#     return FILES


#### MAIN RULE

rule targets:
    input:
        get_genome_meta_files,
        get_genome_meta_files_specialists,
        get_genome_from_ncbi,
        # get_BGCs,
        # get_BGCs_gff,
        get_protein,
        get_amr

# Task 1: Get Generalists metadata
####################################################################################################################

# Import generalists genus
generalists = pd.read_csv("raw/generalists/generalist.txt", header = 0, names = ["generalist"], sep = "\t")
genus = list(generalists["generalist"])
#genus = genus[2:5]

# Cache directories
CACHE_META = "11-generalists-genome-metadata"

# Retrieve the metadata for the generalist genomes
# :q is used to quote
rule genome_metadata_generalists:
    output:
        os.path.join("cache", CACHE_META, "{GENUS}.jsonl")
    conda:
        "src/conda/ncbi-datasets-cli.yaml"
    threads: 1
    shell:
        """
        datasets summary genome taxon {wildcards.GENUS:q} --annotated > {output:q}
        """

# Task 2: Get specialists metadata 
####################################################################################################################

# Import specialist genus
specialists = pd.read_csv("raw/specialists/specialists_mod.txt", header = 0, names = ["specialists"], sep = "\t")
special = list(specialists["specialists"])
#genus = genus[2:5]

# Cache directories
CACHE_SPECIAL = "12-specialists-genome-metadata"

# Retrieve the metadata for the specialists genomes
# :q is used to quote
rule genome_metadata_specialists:
    output:
        os.path.join("cache", CACHE_SPECIAL, "{SPECIAL}.jsonl")
    conda:
        "src/conda/ncbi-datasets-cli.yaml"
    threads: 1
    shell:
        """
        datasets summary genome taxon {wildcards.SPECIAL:q} --annotated > {output:q}
        """

# Task 3: Download genome
####################################################################################################################

# Read tsv file for bacteria
bac_genome = pd.read_csv("raw/generalist_specialist_genome/bacteria.tsv", sep="\t")
b_genome = list(bac_genome["accession"])
#b_genome = b_genome[1:2]


# Read tsv file for bacteria
fun_genome = pd.read_csv("raw/generalist_specialist_genome/fungi.tsv", sep="\t")
f_genome = list(fun_genome["accession"])
# f_genome = f_genome[1:2]

# Perform grouping and then convert to dict
# genomes_dict = genomes.groupby('accession')['kingdom'].apply(list).to_dict()

CACHE_GENOME = "13-generalists-specialists-genome"
CACHE_GFF3 = "13-generalists-specialists-gff3"

rule download_genome:
    output:
        os.path.join("cache", CACHE_GENOME, "{ID}_{KING}_{PRED}.fna"),
        os.path.join("cache", CACHE_GFF3, "{ID}_{KING}_{PRED}.gff")
    conda:
        "src/conda/ncbi-datasets-cli.yaml"
    threads: 1
    params:
        gen_dir = "cache/"+CACHE_GENOME,
        annot_dir = "cache/"+CACHE_GFF3
    shell:
        """
        datasets \
            download genome \
            accession {wildcards.ID} \
            --include genome,gff3 \
            --filename {params.gen_dir}/{wildcards.ID}.zip

        unzip \
            -o {params.gen_dir}/{wildcards.ID}.zip \
            -d {params.gen_dir}
        
        mv \
            {params.gen_dir}/ncbi_dataset/data/{wildcards.ID}/{wildcards.ID}* \
            {params.gen_dir}/{wildcards.ID}_{wildcards.KING}_{wildcards.PRED}.fna
        
        mkdir -p {params.annot_dir}

        mv \
            {params.gen_dir}/ncbi_dataset/data/{wildcards.ID}/*.gff \
            {params.annot_dir}/{wildcards.ID}_{wildcards.KING}_{wildcards.PRED}.gff
        
        rm {params.gen_dir}/{wildcards.ID}.zip
        """

CACHE_BGC = "14-antismash-results"

rule run_antismash_fasta:
    input:
        os.path.join("cache", CACHE_GENOME, "{ID}_{KING}_{PRED}.fna")
    output:
        os.path.join("cache", CACHE_BGC, "{ID}/{ID}_{KING}_{PRED}.json")
    singularity:
        "docker://antismash/standalone:6.1.1"
    threads: 8
    params:
        "cache/"+CACHE_BGC
    shell:
        """
        antismash \
        --taxon {wildcards.KING} \
        --cb-general \
        --cb-knownclusters \
        --cb-subclusters \
        --asf \
        --pfam2go \
        --cpus {threads} \
        --output-dir {params}/{wildcards.ID} \
        --output-basename {wildcards.ID}_{wildcards.KING}_{wildcards.PRED} \
        --genefinding-tool {wildcards.PRED}\
        {input}
        """

CACHE_BGC2 = "14-antismash-results-gff"

rule run_antismash_gff:
    input:
        fa = os.path.join("cache", CACHE_GENOME, "{ID}_{KING}_{PRED}.fna"),
        gff = os.path.join("cache", CACHE_GFF3, "{ID}_{KING}_{PRED}.gff")
    output:
        os.path.join("cache", CACHE_BGC2, "{ID}/{ID}_{KING}_{PRED}.json")
    singularity:
        "docker://antismash/standalone:6.1.1"
    threads: 8
    params:
        "cache/"+CACHE_BGC2
    shell:
        """
        antismash \
        --taxon {wildcards.KING} \
        --cb-general \
        --cb-knownclusters \
        --cb-subclusters \
        --asf \
        --pfam2go \
        --cpus {threads} \
        --output-dir {params}/{wildcards.ID} \
        --output-basename {wildcards.ID}_{wildcards.KING}_{wildcards.PRED} \
        --genefinding-gff3 {input.gff} {input.fa}
        """
# Task 4: Get protein sequence
################################################################################

CACHE_PROTEIN = "13-generalists-specialists-protein"

rule get_proteins:
    output:
        os.path.join("cache", CACHE_PROTEIN, "{ID}_{KING}_{PRED}.fna")
    conda:
        "src/conda/ncbi-datasets-cli.yaml"
    threads: 1
    params:
        gen_dir = "cache/"+CACHE_PROTEIN
    shell:
        """
        datasets \
            download genome \
            accession {wildcards.ID} \
            --include protein \
            --filename {params.gen_dir}/{wildcards.ID}.zip

        unzip \
            -o {params.gen_dir}/{wildcards.ID}.zip \
            -d {params.gen_dir}
        
        mv \
            {params.gen_dir}/ncbi_dataset/data/{wildcards.ID}/protein.faa \
            {params.gen_dir}/{wildcards.ID}_{wildcards.KING}_{wildcards.PRED}.fna
        
        rm {params.gen_dir}/{wildcards.ID}.zip
        """

# Task 5: AMR analysis on bacterial genomes
################################################################################

CACHE_AMR = "16-amr-gff"

rule find_amr:
    input:
        nuc = os.path.join("cache", CACHE_GENOME, "{ID}_bacteria_{PRED}.fna"),
        prot = os.path.join("cache", CACHE_PROTEIN, "{ID}_bacteria_{PRED}.fna"),
        gff = os.path.join("cache", CACHE_GFF3, "{ID}_bacteria_{PRED}.gff")
    output:
        os.path.join("cache", CACHE_AMR, "{ID}_bacteria_{PRED}_amr.tsv")
    singularity:
        "docker://ncbi/amr:3.11.8-2023-02-23.1"
    threads: 8
    shell:
        """
        amrfinder \
            --nucleotide {input.nuc} \
            --protein {input.prot} \
            --gff {input.gff} \
            --plus \
            --threads {threads} \
            -o {output}
        """



