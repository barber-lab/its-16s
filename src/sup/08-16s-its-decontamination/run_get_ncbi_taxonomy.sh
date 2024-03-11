#!/bin/bash

# Activate local conda environment
source /home/${USER}/.bashrc
module purge
mamba activate python3.8

WD="/home/qi47rin/proj/04-global-microbiome"
BASE="/home/qi47rin/proj/02-compost-microbes"

cd $WD

mkdir -p cache/08-contamination/03-taxonomy-from-ncbi

python $BASE/src/07-fungal-bac-abundances/get_taxonomy_from_ncbi.py \
    --input cache/08-contamination/02-taxonomy-id-list/contaminants.tsv \
    --output cache/08-contamination/03-taxonomy-from-ncbi/tax.tsv
