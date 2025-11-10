#!/usr/bin/env bash
#SBATCH -J bam_coverage_analysis
#SBATCH --ntasks=10
#SBATCH -N 1
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem=100G
#SBATCH -o /scratch/ejy4bu/erroroutputs/bam_analysis/coverage.%A.out
#SBATCH -e /scratch/ejy4bu/erroroutputs/bam_analysis/coverage.%A.err
#SBATCH -p standard
#SBATCH --account berglandlab

mkdir -p /scratch/ejy4bu/erroroutputs/bam_analysis

#module load samtools
module load bedtools

#https://www.htslib.org/doc/samtools-depth.html
# https://www.htslib.org/doc/samtools-coverage.html

# make windows once:
# head -n 14 /scratch/ejy4bu/compBio/genomefiles/scaffold_lengths.txt > /scratch/ejy4bu/compBio/genomefiles/chr1-14_names_and_lengths.txt
# chr_list="/scratch/ejy4bu/compBio/genomefiles/chr1-14_names_and_lengths.txt"
# bedtools makewindows -g "$chr_list" -w 500 > /scratch/ejy4bu/compBio/genomefiles/windows_500bp.bed
windows_500="/scratch/ejy4bu/compBio/genomefiles/windows_500bp.bed"

bam_dir="/scratch/ejy4bu/compBio/bams"
out_dir="/scratch/ejy4bu/compBio/bam_analysis/coverage"
mkdir -p "$out_dir"

#bam="/scratch/ejy4bu/compBio/bams/Robert_samples_bams/RobertUK_F1/RobertUK_F1.dedup.bam"

for bam in "$bam_dir"/*/*/*.dedup.bam; do
    sample=$(basename "$bam" .dedup.bam)
    mkdir -p "$out_dir/${sample}"
    out_file="$out_dir/${sample}/${sample}_500bp.tsv"

    echo -e "chrom\tstart\tend\tmean_depth" > "$out_file"
    echo "Processing $sample..."

    bedtools coverage -a "$windows_500" -b "$bam" -mean >> "$out_file"

    tr '\t' ',' < "$out_file" > "$out_dir/${sample}/${sample}_500bp.csv"
done


### generates read for each base pair

# echo -e "chrom\pos\depth" > "$out_file"

# while read -r chr; do
#     echo "Processing $chr..."
#     samtools depth -r "$chr" "$bam" >> "$out_file"
# done < $chr_list

# tr '\t' ',' < /scratch/ejy4bu/compBio/bam_analysis/chr1_depth.tsv > /scratch/ejy4bu/compBio/bam_analysis/chr1_depth.csv
# rm -f $outfile


