#!/usr/bin/env bash
#
#SBATCH -J run_kraken # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/erroroutputs/down.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/erroroutputs/down.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

#ijob -A berglandlab -c10 -p standard --mem=40G

module load kraken2

DBNAME="/scratch/ejy4bu/compBio/kraken"

#standard kraken db 
# kraken2-build --standard --threads 24 --db $DBNAME

#custom database
# kraken2-build --download-taxonomy --db $DBNAME
kraken2-build --threads 10 --download-library nt --db $DBNAME

:<<comment

kraken2-build --build --threads 10 --db $DBNAME



#robert's method:
kraken2-build --build --threads 10 --download-library nt --db nt


kraken2-build --build --threads 10 -db --download-library nt 

kraken2-build --standard --threads 10 --db nt 
kraken2 --db $DBNAME seqs.fa

kraken2 --db $DBNAME seqs.fa

# do once

# download ft library
# unzip it
# build 
DBNAME="/scratch/ejy4bu/compBio/kraken"
kraken2-build --threads 10 --download-library nt --db $DBNAME
tar 
kraken2-build --standard --threads 10 --db $DBNAME

REPORTS="/scratch/ejy4bu/compBio/kraken/reports

kraken2 --db $DBNAME \
    --threads 10 --report $REPORTS\  \
    --classified-out 
    --output
    --use-names
    $sample

comment