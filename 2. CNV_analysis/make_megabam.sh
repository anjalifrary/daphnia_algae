#!/usr/bin/env bash
#
#SBATCH -J megabam 
#SBATCH --array=1-3
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/err_outs/megabam/megabam.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/megabam/megabam.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load samtools

out_dir="/scratch/ejy4bu/compBio/cnv/megabams"
mkdir -p $out_dir

groups=(/scratch/ejy4bu/compBio/cnv/data_tables/bam_paths/*_bamList.txt)
bam_list=${groups[$SLURM_ARRAY_TASK_ID-1]}
group=$(basename $bam_list _bamList.txt)
#in_bams="/scratch/ejy4bu/compBio/cnv/data_tables/bam_paths/${group}_bamList.txt"
out_bam="${out_dir}/${group}.megabam.bam"

samtools merge -f -@ 10 $out_bam -b $bam_list
samtools index $out_bam 

echo "Done merging $group"


# megabam="/scratch/ejy4bu/compBio/cnv/megabams/REED_NotSephedex.megabam.bam"
# samtools quickcheck -v $megabam
# samtools view -H $megabam | head
# samtools idxstats $megabam
# samtools flagstat $megabam
