#!/usr/bin/env bash
#SBATCH -J cn.mops
#SBATCH --ntasks=10
#SBATCH -N 1
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem=100G
#SBATCH -o /scratch/ejy4bu/err_outs/cnmops/cnmops.%A.out
#SBATCH -e /scratch/ejy4bu/err_outs/cnmops/cnmops.%A.err
#SBATCH -p standard
#SBATCH --account berglandlab

mkdir -p /scratch/ejy4bu/err_outs/cnmops

export R_LIBS_USER="/sfs/gpfs/tardis/home/ejy4bu/R/goolf/4.5/"
module load gcc/11.4.0 openmpi/4.1.4 icu R/4.5.0

module load R/4.5.0

Rscript 2.\ CNV_analysis/cn-mops.R