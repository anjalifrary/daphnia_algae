#!/usr/bin/env bash
#
#SBATCH -J megabam 
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/err_outs/megabam/megabam.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/megabam/megabam.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

out_dir="/scratch/ejy4bu/compBio/cnv/megabams"
mkdir -p $out_dir

group="REED_NotSephedex"
in_bams="/scratch/ejy4bu/compBio/cnv/data_tables/bam_paths/${group}_bamList.txt"
out_bam="${out_dir}/${group}.megabam.bam"

samtools merge -f -@ 10 $out_bam -b $in_bams 
samtools index $out_bam 
