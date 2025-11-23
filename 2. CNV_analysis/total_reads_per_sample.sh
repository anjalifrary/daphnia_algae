#!/usr/bin/env bash
#
#SBATCH -J subset 
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/err_outs/subset/subset.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/subset/subset.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load samtools

bam_dir="/scratch/ejy4bu/compBio/cnv/bams"
out_file="/scratch/ejy4bu/compBio/cnv/bam_read_counts.csv"

# write header to csv output file
echo "sample,group,total_reads" > $out_file

for bam in ${bam_dir}/*/*/*.dedup.bam; do
    sample=$(basename $bam .dedup.bam)
    group=$(basename $(dirname $(dirname $bam)))

    echo "Processing sample $sample in group $group"

    # count total mapped reads, -F 4 exclude unmapped reads
    reads=$(samtools view -c -F 4 "$bam")
    echo "${sample},${group},${reads}" >> $out_file
done

echo "Finished calculating reads for all samples"