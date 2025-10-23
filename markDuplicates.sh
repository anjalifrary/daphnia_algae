#!/bin/bash
#SBATCH --job-name=markdups         # Job name
#SBATCH --output=/scratch/ejy4bu/erroroutputs/markdups_%j.out
#SBATCH --error=/scratch/ejy4bu/erroroutputs/markdups_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=04:00:00
#SBATCH -p standard
#SBATCH --account berglandlab
#SBATCH --array=1-320%20

module load picard

bams=(/scratch/ejy4bu/compBio/bams/*/*/*.sort.bam)

bam=${bams[$((SLURM_ARRAY_TASK_ID - 1))]}
sample=$(basename "$bam" .sort.bam)
outdir=$(dirname "$bam")
dedup_bam="${outdir}/${sample}.dedup.bam"


echo "Processing $bam"

if [ -f "$dedup_bam" ]; then
    echo "Skipping $bam (already deduplicated)"
    exit 0
fi

java -Xmx45G -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    I="$bam" \
    O="$dedup_bam" \
    M="${outdir}/${sample}.dedup.metrics" \
    REMOVE_DUPLICATES=true \
    CREATE_INDEX=true

# for bam in /scratch/ejy4bu/compBio/bams/*/*/*.sort.bam; do
#     sample=$(basename $bam .sort.bam)
#     outdir=$(dirname $bam)
#     java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
#         I=$bam \
#         O=${outdir}/${sample}.dedup.bam \
#         M=${outdir}/${sample}.dedup.metrics \
#         REMOVE_DUPLICATES=true \
#         CREATE_INDEX=true
# done