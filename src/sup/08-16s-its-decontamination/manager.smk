import pandas as pd
import os
configfile: "/home/qi47rin/proj/04-global-microbiome/src/08-decontamination/smk.config"
workdir: config['workdir']

module get_amplicon_pe:
    snakefile: "/home/qi47rin/proj/04-global-microbiome/src/00-pipelines/02-pe-fq-gz-sra.smk"
    config: config

module get_amplicon_se:
    snakefile: "/home/qi47rin/proj/04-global-microbiome/src/00-pipelines/02-se-fq-gz-sra.smk"
    config: config

module quality_control_pe:
    snakefile: "/home/qi47rin/proj/04-global-microbiome/src/00-pipelines/05-quality_control_pe.smk"
    config: config

module quality_control_se:
    snakefile: "/home/qi47rin/proj/04-global-microbiome/src/00-pipelines/05-quality_control_se.smk"
    config: config

module taxonomic_profiling:
    snakefile: "/home/qi47rin/proj/04-global-microbiome/src/00-pipelines/01-closed-otu-picking.smk"
    config: config

# When running the function, the IDS become available
def compressed_raw_fq_pe(file_path):
    pe_tbl = pd.read_csv(file_path, header = 0, names = ["Run"], sep = "\t")
    global pe_ids
    pe_ids = list(pe_tbl["Run"])
    PE = expand(os.path.join(config['PROJ'], config['FOLDERS']['raw_fq_pe'], "{ID}_{PE}.fq.gz"), ID = pe_ids, PE = [1,2])
    return PE

def compressed_raw_fq_se(file_path):
    se_tbl = pd.read_csv(file_path, header = 0, names = ["Run"], sep = "\t")
    global se_ids
    se_ids = list(se_tbl["Run"])
    SE = expand(os.path.join(config['PROJ'], config['FOLDERS']['raw_fq_se'], "{ID}.fq.gz"), ID = se_ids)
    return SE

use rule * from get_amplicon_pe

use rule * from get_amplicon_se

use rule * from quality_control_pe

use rule * from quality_control_se

use rule * from taxonomic_profiling

rule all:
    input:
        compressed_raw_fq_pe(config['SAMPLES_PE']),
        compressed_raw_fq_se(config['SAMPLES_SE']),
        os.path.join(config['PROJ'], config['FOLDERS']['stats']),
        expand(os.path.join(config['PROJ'], config['FOLDERS']['manifest'], "{DB}.tsv"), DB=config['PARAMETERS']['database']),
        expand(os.path.join(config['PROJ'], config['FOLDERS']['fastqc_trim_merg'], "{ID}_fastqc.html"), ID = pe_ids),
        expand(os.path.join(config['PROJ'], config['FOLDERS']['fastqc_trim_se'], "{ID}_fastqc.html"), ID = se_ids),
        expand(os.path.join(config['PROJ'], config['FOLDERS']['reference_artifact'], "{DB}_seq.qza"), DB=config['PARAMETERS']['database']),
        expand(os.path.join(config['PROJ'], config['FOLDERS']['abundance'], "{DB}-feature-table-abundances.tsv"), DB=config['PARAMETERS']['database'])
    default_target: True


rule multiqc_decontam:
    input:
        raw_pe = expand(os.path.join(config['PROJ'], config['FOLDERS']['fastqc_raw_pe'], "{ID}_{PE}_fastqc.html"), ID = pe_ids, PE = [1,2]),
        raw_se = expand(os.path.join(config['PROJ'], config['FOLDERS']['fastqc_raw_se'], "{ID}_fastqc.html"), ID = se_ids),
        trim_merg = expand(os.path.join(config['PROJ'], config['FOLDERS']['fastqc_trim_merg'], "{ID}_fastqc.html"), ID = pe_ids),
        se_trim = expand(os.path.join(config['PROJ'], config['FOLDERS']['fastqc_trim_se'], "{ID}.html"), ID = se_ids)
    output:
        directory(os.path.join(config['PROJ'], config['FOLDERS']['multiqc']))
    params:
        raw_pe = os.path.join(config['PROJ'], config['FOLDERS']['fastqc_raw_pe']),
        raw_se = os.path.join(config['PROJ'], config['FOLDERS']['fastqc_raw_se']),
        trim_merg = os.path.join(config['PROJ'], config['FOLDERS']['fastqc_trim_merg']),
        se_trim = os.path.join(config['PROJ'], config['FOLDERS']['fastqc_trim_se'])
    singularity:
        config['IMAGES']['multiqc']
    shell:
        """
        multiqc \
            --force \
            --dirs \
            --dirs-depth 1 \
            --outdir {output} {params.raw_pe} {params.raw_se} {params.se_trim} {params.trim_merg}
        """

rule classify_runs:
    input:
        link=compressed_raw_fq_pe(config['SAMPLES_PE']),
        script=config['SCRIPTS']['classify']
    output:
        directory(os.path.join(config['PROJ'], config['FOLDERS']['stats']))
    params:
        "src/08-decontamination"
    shell:
        """
        bash {input.script}
        """

rule silva_unite_manifest:
    input:
        link=rules.classify_runs.output,
        script=config['SCRIPTS']['manifest']
    output:
        os.path.join(config['PROJ'], config['FOLDERS']['manifest'], "{DB}.tsv"),
    conda:
        os.path.join(config['workdir'], config['ENVIRONMENTS']['rtidyverse'])
    resources:
        mem_mb=100000
    shell:
        """
        Rscript {input.script}
        """
