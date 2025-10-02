#!/usr/bin/env bash
#
#SBATCH -J run_kraken # A single job name for the array
#SBATCH --ntasks-per-node=4 # ten core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 300G
#SBATCH -o /scratch/ejy4bu/erroroutputs/down.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/erroroutputs/down.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


module load kraken2

DBNAME="/scratch/ejy4bu/compBio/kraken/nt"

export KRAKEN2_DATA_PATH="/scratch/ejy4bu/compBio/kraken/nt"

#kraken2-build --db $DBNAME --build --threads 10

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
samp_path="/scratch/ejy4bu/compBio/Robert_samples/RobertUK_G12"
sample="RobertUK_G12_CKDL250003065-1A_22M5YKLT4_L4"

SAMPLE1="${samp_path}/${sample}_1.P.trimm.fastq"
SAMPLE2="${samp_path}/${sample}_2.P.trimm.fastq"

# take 10k random
#seqtk sample -s 100 "$SAMPLE1" 1000 > ${samp_path}/${sample}_1.P.trimm.fastq
#seqtk sample -s 100 "$SAMPLE2" 1000 > ${samp_path}/${sample}_2.P.trimm.fastq

echo "Running samples ${sample}"

:<<sample_def
#long read:
SAMPLE="/project/berglandlab/chlorella_sequencing/raw_longread_from_Reed/m84128_250121_222443_s2.hifi_reads.bc2104.fq.gz"

#reference:
SAMPLE="/project/berglandlab/chlorella_sequencing/reference_genome/GCA_023343905.1_cvul_genomic.fa"
sample_def

#report folder
REPORTS="/scratch/ejy4bu/compBio/kraken/reports"
mkdir -p "${REPORTS}" 

#:<<paired
kraken2 --db $DBNAME \
    --threads 4 \
    --output ${REPORTS}/${sample}_output.txt \
    --report ${REPORTS}/${sample}_report.txt \
    --classified-out ${REPORTS}/${sample}_#_classified.fq \
    --use-names \
    --paired "${samp_path}/${sample}/${sample}_1_sub.fastq" "${samp_path}/${sample}/${sample}_2_sub.fastq"

#paired

:<<unpaired
kraken2 --db $DBNAME \
    --threads 4 \
    --output ${REPORTS}/$(basename "${SAMPLE}")_output.txt \
    --report ${REPORTS}/$(basename "${SAMPLE}")_report.txt  \
    --classified-out ${REPORTS}/$(basename "${SAMPLE}")_classified.fq \
    --use-names \
    $SAMPLE
unpaired