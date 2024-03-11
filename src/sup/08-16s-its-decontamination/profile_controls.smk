configfile: "/home/qi47rin/proj/04-global-microbiome/src/08-decontamination/smk.config"

import pandas as pd
import os
workdir: config['workdir']


module taxonomic_profiling:
    snakefile: "/home/qi47rin/proj/04-global-microbiome/src/00-pipelines/01-closed-otu-picking.smk"
    config: config

use rule * from taxonomic_profiling

# When running the function, the IDS become available
def compressed_raw_fq_pe(file_path):
    pe_tbl = pd.read_csv(file_path, header = 0, names = ["Run"], sep = "\t")
    global pe_ids
    pe_ids = list(pe_tbl["Run"])
    PE = expand(os.path.join(config['PROJ'], config['FOLDERS']['raw_fq_pe'], "{ID}_{PE}.fq.gz"), ID = pe_ids, PE = [1,2])
    return PE

use rule * from alignment

use rule * from count_aligned

rule all:
    input:
        compressed_raw_fq_pe(config['REFERENCE']['fq_pe']),
        os.path.join(config['PROJ'], config['GENOME_INDEX']),
        expand(os.path.join(config['PROJ'], config['FOLDERS']['sorted'], "{ID}.sorted.bam"), ID = pe_ids),
        expand(os.path.join(config['PROJ'], config['FOLDERS']['stats'], "{ID}.txt"), ID = pe_ids)
    default_target: True

