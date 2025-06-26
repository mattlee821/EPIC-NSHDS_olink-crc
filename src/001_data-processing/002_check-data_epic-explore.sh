#!/bin/bash

#SBATCH --job-name=002_check-data_epic-explore
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00-03:00:00
#SBATCH --mem=50000M
#SBATCH --partition=low_p

export TMPDIR=/scratch/leem/temp/

cd /data/IARC_Biostat/work/EPIC-NSHDS_olink-crc/

VAR1=002_check-data_epic-explore
/opt/R/4.4.1/bin/Rscript src/001_data-processing/${VAR1}.R
