#!/bin/bash -l
#SBATCH -A grp-org-gt-seqtech
#SBATCH -p dori
#SBATCH -q jgi_long
#SBATCH -t 200:00:00
#SBATCH -c 4
#SBATCH --mem=50G
#SBATCH --job-name=plot_epimark
#SBATCH -o plot_epimark.log
#SBATCH -e plot_epimark.err

cd /clusterfs/jgi/groups/gentech/homes/romalley/full_DAPseq_annotation/script

source activate r-kernel
# R CMD BATCH 10_plot_epimark_distribution_TF.R
R CMD BATCH 10_plot_epimark_distribution_TF_intersect.R
