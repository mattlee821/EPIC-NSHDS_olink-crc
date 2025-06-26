#!/bin/bash

#SBATCH --job-name=olink-heterogeneity
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00-10:00:00
#SBATCH --mem=80000M
#SBATCH --partition=low_p

VAR1=001_master-parallel

export TMPDIR=/scratch/leem/temp/ 

cd /data/IARC_Biostat/work/EPIC-NSHDS_olink-crc/

/opt/R/4.4.1/bin/Rscript src/002_EPIC-analysis/002_heterogeneity/${VAR1}.R
