#!/bin/bash

echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $(pwd)"
echo ""
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"
echo "" 

# Activate local conda environment
source /home/${USER}/.bashrc

SCRATCH="/scratch/qi47rin"

mamba activate snakemake

cd src/08-decontamination

snakemake  \
    --snakefile classify_samples.smk \
    --profile /home/qi47rin/proj/04-global-microbiome/src/ \
    --singularity-prefix cache/00-singularity \
    --conda-prefix cache/00-conda-env \
    --singularity-args "--bind $SCRATCH" \
    --conda-frontend mamba \
    --nolock
