#!/usr/bin/env bash
#
#SBATCH -J freec
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00 # 1 day
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/err_outs/freec/freec.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/freec/freec.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# make conda environment with freec in it

#conda create -n freec-env -c bioconda control-freec python=3.10
conda activate freec-env

### failed because conda version of freec is old
#   and download links don't include necessary libraries


