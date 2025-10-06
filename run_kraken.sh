#!/usr/bin/env bash
#
#SBATCH -J run_kraken # A single job name for the array
#SBATCH --ntasks-per-node=4 # ten core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 300G
#SBATCH -o /scratch/ejy4bu/erroroutputs/kraken.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/erroroutputs/kraken.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab


module load kraken2
module load samtools

DBNAME="/scratch/ejy4bu/compBio/kraken/nt"

export KRAKEN2_DATA_PATH="/scratch/ejy4bu/compBio/kraken/nt"

#kraken2-build --db $DBNAME --build --threads 10

# kraken2-build --download-library nt --threads 10 --db $DBNAME

#standard kraken db 
# kraken2-build --standard --threads 24 --db $DBNAME

#custom database
# kraken2-build --download-taxonomy --db $DBNAME

samp_path="/scratch/ejy4bu/compBio/Robert_samples_bams"
sample="$1"
#sample="${samp_path}/RobertUK_F1.sort.bam"
samp_name=$(basename ${sample%.sort.bam})

SAMPLE="${samp_path}/${samp_name}/${samp_name}_output.fastq"


#convert bam back to fastq
#samtools fastq -o ${SAMPLE} ${sample}
if [ ! -f "$SAMPLE" ]; then
    echo "No fastq file found"
    exit 1
fi

:<<sample_def
#short read:
samp_path="/scratch/ejy4bu/compBio/Robert_samples/RobertUK_G12"
sample="RobertUK_G12_CKDL250003065-1A_22M5YKLT4_L4"

SAMPLE1="${samp_path}/${sample}_1.P.trimm.fastq"
SAMPLE2="${samp_path}/${sample}_2.P.trimm.fastq"

# take 10k random
# seqtk sample -s 100 "$SAMPLE1" 1000 > ${samp_path}/${sample}_1.P.trimm.sub.fastq
# seqtk sample -s 100 "$SAMPLE2" 1000 > ${samp_path}/${sample}_2.P.trimm.sub.fastq

echo "Running samples ${sample}"


#long read:
sample="/project/berglandlab/chlorella_sequencing/raw_longread_from_Reed/m84128_250121_222443_s2.hifi_reads.bc2104.fq.gz"

#reference:
sample="/project/berglandlab/chlorella_sequencing/reference_genome/GCA_023343905.1_cvul_genomic.fa"
sample_def

#report folder
REPORTS="/scratch/ejy4bu/compBio/kraken/reports/${samp_name}"
mkdir -p "${REPORTS}" 

echo "Running kraken on ${samp_name}"

:<<paired
kraken2 --db $DBNAME \
    --threads 4 \
    --output ${REPORTS}/${sample}_output.txt \
    --report ${REPORTS}/${sample}_report.txt \
    --classified-out ${REPORTS}/${sample}_#_classified.fq \
    --use-names \
    --paired "${samp_path}/${sample}_1.P.trimm.sub.fastq" "${samp_path}/${sample}_2.P.trimm.sub.fastq"

paired

#:<<unpaired
kraken2 --db $DBNAME \
    --threads 4 \
    --output ${REPORTS}/$(basename "${SAMPLE}")_output.txt \
    --report ${REPORTS}/$(basename "${SAMPLE}")_report.txt  \
    --classified-out ${REPORTS}/$(basename "${SAMPLE}")_classified.fq \
    --use-names \
    $SAMPLE
#unpaired