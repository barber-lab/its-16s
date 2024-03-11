configfile: "/home/qi47rin/proj/04-global-microbiome/src/08-decontamination/smk.config"

import pandas as pd
import os
workdir: config['workdir']

module alignment_pe:
    snakefile: "/home/qi47rin/proj/04-global-microbiome/src/00-pipelines/03-alignment_pe.smk"
    config: config

module alignment_se:
    snakefile: "/home/qi47rin/proj/04-global-microbiome/src/00-pipelines/03-alignment_se.smk"
    config: config

module count_aligned:
    snakefile: "/home/qi47rin/proj/04-global-microbiome/src/00-pipelines/04-counts-from-alignment-NW.smk"
    config: config

# When running the function, the IDS become available
def sorted_bam_files_pe(file_path):
    pe_tbl = pd.read_csv(file_path, header = 0, names = ["Run"], sep = "\t")
    global pe_ids
    pe_ids = list(pe_tbl["Run"])
    PE = expand(os.path.join(config['PROJ'], config['FOLDERS']['sorted'], "{ID}_{BASE_NAME}.sorted.bam"), ID = pe_ids, BASE_NAME = config['PARAMETERS']['database'])
    return PE

def sorted_bam_files_se(file_path):
    se_tbl = pd.read_csv(file_path, header = 0, names = ["Run"], sep = "\t")
    global se_ids
    se_ids = list(se_tbl["Run"])
    SE = expand(os.path.join(config['PROJ'], config['FOLDERS']['sorted'], "{ID}_{BASE_NAME}.sorted.bam"), ID = se_ids, BASE_NAME = config['PARAMETERS']['database'])
    return SE

use rule * from alignment_pe

use rule * from alignment_se

use rule * from count_aligned

rule all:
    input:
        expand(os.path.join(config['REFERENCE'], "{BASE_NAME}.fna"), BASE_NAME = config['PARAMETERS']['database']),
        sorted_bam_files_pe(config['SAMPLES_PE']),
        sorted_bam_files_se(config['SAMPLES_SE']),
        expand(os.path.join(config['PROJ'], config['FOLDERS']['stats'], "{ID}_{BASE_NAME}.txt"), ID = pe_ids+se_ids, BASE_NAME = config['PARAMETERS']['database']),
    default_target: True

