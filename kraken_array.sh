#!/usr/bin/env bash
#SBATCH -J kraken_array # A single job name for the array
#SBATCH --array=1-25
#SBATCH --ntasks-per-node=4 # 4 core
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH -N 1 # on one node
#SBATCH -t 0-12:00 # 10 hours
#SBATCH --mem 300G
#SBATCH -o /scratch/ejy4bu/erroroutputs/kraken.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/erroroutputs/kraken.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

BAM_DIR="/scratch/ejy4bu/compBio/Robert_samples_bams"

mapfile -t BAM_FILES < <(find "$BAM_DIR" -mindepth 2 -maxdepth 2 -name "*.sort.bam" | sort)

BAM="${BAM_FILES[$SLURM_ARRAY_TASK_ID-1]}"

samp_name=$(basename ${BAM%.sort.bam})
FASTQ="${BAM_DIR}/${samp_name}/${samp_name}_output.fastq"

if [ ! -f "$FASTQ" ]; then
    echo "Converting bam $BAM to fastq $FASTQ"
    samtools fastq -@ 10 -o ${FASTQ} ${BAM}
else 
    echo "Fastq already exists $FASTQ"
fi

echo "Array task $SLURM_ARRAY_TASK_ID processing BAM: $BAM"

run_kraken.sh "$BAM"
