#!/usr/bin/env bash
#SBATCH -J run_copy_number_variants
#SBATCH --ntasks=10
#SBATCH -N 1
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem=100G
#SBATCH -o /scratch/ejy4bu/erroroutputs/bam_analysis/cnv.%A.out
#SBATCH -e /scratch/ejy4bu/erroroutputs/bam_analysis/cnv.%A.err
#SBATCH -p standard
#SBATCH --account berglandlab

mkdir -p /scratch/ejy4bu/erroroutputs/bam_analysis

export R_LIBS_USER="/sfs/gpfs/tardis/home/ejy4bu/R/goolf/4.5/"
module load gcc/11.4.0  openmpi/4.1.4 icu R/4.5.0

Rscript /home/ejy4bu/daphnia_algae/copy_number_variants_analysis.R

