#!/usr/bin/env bash
#
#SBATCH -J kraken_array # A single job name for the array
#SBATCH --ntasks-per-node=4 # 4 core
#SBATCH -N 1 # on one node
#SBATCH -t 0-12:00 # 10 hours
#SBATCH --mem 300G
#SBATCH -o /scratch/ejy4bu/erroroutputs/kraken.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/erroroutputs/kraken.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-