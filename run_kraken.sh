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


module load kraken2

DBNAME="/scratch/ejy4bu/compBio/kraken/nt"

kraken2-build --db $DBNAME --build --threads 10

# kraken2-build --download-library nt --threads 10 --db $DBNAME

#standard kraken db 
# kraken2-build --standard --threads 24 --db $DBNAME

#custom database
# kraken2-build --download-taxonomy --db $DBNAME

:<<comment

kraken2-build --build --threads 10 --db $DBNAME

#robert's method:
kraken2-build --build --threads 10 --download-library nt --db nt


kraken2-build --build --threads 10 -db --download-library nt 

kraken2-build --standard --threads 10 --db nt 
kraken2 --db $DBNAME seqs.fa

kraken2 --db $DBNAME seqs.fa

# do once

# download nt library
# unzip it
# build 
DBNAME="/scratch/ejy4bu/compBio/kraken"
kraken2-build --threads 10 --download-library nt --db $DBNAME
tar 
kraken2-build --standard --threads 10 --db $DBNAME
comment


#short read:
SAMPLE1="/scratch/ejy4bu/compBio/fastq/SRR14426881_1.fastq"
SAMPLE2="/scratch/ejy4bu/compBio/fastq/SRR14426881_2.fastq"
SAMPLE=$(basename "${SAMPLE1}" _1.fastq)

echo "Running samples ${SAMPLE}"

:<<sample_def
#long read:
SAMPLE="/project/berglandlab/chlorella_sequencing/raw_longread_from_Reed/m84128_250121_222443_s2.hifi_reads.bc2104.fq.gz"

#reference:
SAMPLE="/project/berglandlab/chlorella_sequencing/reference_genome/GCA_023343905.1_cvul_genomic.fa"
sample_def

#report folder
REPORTS="/scratch/ejy4bu/compBio/kraken/reports"
mkdir -p "${REPORTS}" 

kraken2 --db $DBNAME \
    --threads 10 \
    --paired "${SAMPLE1}" "${SAMPLE2}" \
    --output ${REPORTS}/${SAMPLE}_output.txt \
    --report ${REPORTS}/${SAMPLE}_report.txt \
    --classified-out ${REPORTS}/${SAMPLE}_#_classified.fq \
    --use-names \
    --memory-mapping


:<<unpaired
kraken2 --db $DBNAME \
    --threads 10 \
    --fastq-input $SAMPLE \
    --output ${REPORTS}/$(basename "${SAMPLE}")_output.txt \
    --report ${REPORTS}/$(basename "${SAMPLE}")_report.txt  \
    --classified-out ${REPORTS}/$(basename "${SAMPLE}")_classified.fq \
    --use-names
unpaired