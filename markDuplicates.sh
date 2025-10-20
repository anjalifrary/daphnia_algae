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

module load picard

for bam in /scratch/ejy4bu/compBio/bams/*/*/*.sort.bam; do
    sample=$(basename $bam .sort.bam)
    outdir=$(dirname $bam)
    picard MarkDuplicates \
        I=$bam \
        O=${outdir}/${sample}.dedup.bam \
        M=${outdir}/${sample}.dedup.metrics \
        REMOVE_DUPLICATES=true \
        CREATE_INDEX=true
done