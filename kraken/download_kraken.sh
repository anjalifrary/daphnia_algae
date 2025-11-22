#!/usr/bin/env bash
#SBATCH -J download_kraken # A single job name for the array
#SBATCH --ntasks-per-node=2 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 10:00:00 ### 15 seconds
#SBATCH --mem 20G
#SBATCH -o /scratch/ejy4bu/erroroutputs/move.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/erroroutputs/move.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --mail-type=END               # Send email at job completion
#SBATCH --mail-user=ejy4bu@virginia.edu    # Email address for notifications

# cp -r /scratch/ejy4bu/compBio/kraken /project/berglandlab/chlorella_sequencing/krakendbtmp

PROJECT_DB="/project/berglandlab/chlorella_sequencing/krakendbtmp"
MYSCRATCH="/scratch/ejy4bu/compBio/kraken"
mkdir -p $MYSCRATCH
rsync -avh --progress $PROJECT_DB/ $MYSCRATCH/

rm -rf /project/berglandlab/chlorella_sequencing/krakendbtm
