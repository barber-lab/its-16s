#!/bin/bash

#SBATCH --job-name=DECONTAMINATION                                          # Job name
#SBATCH --ntasks=1                                                          # Run a single task
#SBATCH --cpus-per-task=1                                                   # Number of CPU cores per task
#SBATCH --mem=20G                                                           # Job memory request
#SBATCH --time=3-00:00:00                                                   # Time limit hrs:min:sec
#SBATCH --output=../../logs/09-decontamination/decontamination_%j.log       # Standard output and error log
#SBATCH --partition=standard

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

snakemake  \
    --snakefile manager.smk \
    --profile /home/qi47rin/proj/02-compost-microbes/src/ \
    --singularity-prefix cache/00-singularity \
    --conda-prefix cache/00-conda-env \
    --singularity-args "--bind $SCRATCH" \
    --conda-frontend mamba \
    --nolock
