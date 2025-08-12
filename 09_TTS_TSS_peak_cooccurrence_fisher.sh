#!/bin/bash -l
#SBATCH -A grp-org-gt-seqtech
#SBATCH -p dori
#SBATCH -q jgi_long
#SBATCH -t 200:00:00
#SBATCH -c 4
#SBATCH --mem=10G
#SBATCH --job-name=tomtom_parallel
#SBATCH -o TTS_TSS_cooccur.log
#SBATCH -e TTS_TSS_cooccur.err

cd /clusterfs/jgi/groups/gentech/homes/romalley/full_DAPseq_annotation/script

source activate r-kernel
R CMD BATCH 09_TTS_TSS_peak_cooccurrence_fisher.R