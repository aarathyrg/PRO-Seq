#!/usr/bin/env bash

# proseq-2_0 index danko_lab
# Author: Aarathy and Katrin
# Date: 2022-02-14

# got files with

# need to be extract

# === SLURM directives ===

# Run this on one node and as a single task (these are the default settings anyway)
#SBATCH --job-name=proseq2.0
#SBATCH --nodes=1
#SBATCH --array=0-161
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=20G
# Max. runtime is 120 minutes
#SBATCH --time=47:00:00
# The "short" Quality of Time (QoS) is the default
# Set this to "medium" if you need a walltime between 8 hours and 2 days
# and to "long" for jobs taking more than 2 days (but less than 14)
#SBATCH --qos=medium

# Give a sensible name to the output file
# (stdout and stderr both go to this)
#SBATCH --output=proseq2_0_out.txt

FILES=($(ls reads | cut -d "_" -f1,2,3,4))
# === DATA SETUP ===

# Create a working directory under your own scratch directory
WORKDIR=/scratch-cbe/users/aarathy.geetha/PRO-seq_round2_novaseq/

cd $WORKDIR         # Change to the working directory

# LOAD MODULES
ml restore pro_env

# === DO THE WORK ===

bash proseq2.0.bsh -PE -i bwa_index/Mus_musculus.GRCm38.dna.toplevel.fa.gz -c genome.fa.sizes -I reads/${FILES[$SLURM_ARRAY_TASK_ID]} -T tmp -O results  --RNA3=R2_5prime -3 -s  --ADAPT1=TGGAATTCTCGGGTGCCAAGG --ADAPT2=GATCGTCGGACTGTAGAACTCTGAAC --UMI1=4 --UMI2=8 --thread=8


# Congratulations! This is the end of the script!
