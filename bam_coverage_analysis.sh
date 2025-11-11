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
# bedtools makewindows -g "$chr_list" -w 5000 > /scratch/ejy4bu/compBio/genomefiles/windows_5000bp.bed
bedtools makewindows -g "$chr_list" -w 10000 > /scratch/ejy4bu/compBio/genomefiles/windows_10000bp.bed
windows="/scratch/ejy4bu/compBio/genomefiles/windows_10000bp.bed"

bam_dir="/scratch/ejy4bu/compBio/bams"
out_dir="/scratch/ejy4bu/compBio/bam_analysis/coverage_data"
mkdir -p "$out_dir"

#bam="/scratch/ejy4bu/compBio/bams/Robert_samples_bams/RobertUK_F1/RobertUK_F1.dedup.bam"
#bam="/scratch/ejy4bu/compBio/bams/Old_Algae_bams/SRR14370481/SRR14370481.dedup.bam"

for bam in "$bam_dir"/*/*/*.dedup.bam; do
    group=$(basename "${bam%/*/*}")
    sample=$(basename "$bam" .dedup.bam)
    mkdir -p "$out_dir/${group}/${sample}"
    out_file="$out_dir/${group}/${sample}/${sample}_10000bp.tsv"

    echo -e "chrom\tstart\tend\tmean_depth" > "$out_file"
    echo "Processing $sample..."

    bedtools coverage -a "$windows" -b "$bam" -mean >> "$out_file"

    tr '\t' ',' < "$out_file" > "$out_dir/${group}/${sample}/${sample}_10000bp.csv"
done

# Run copy_number_variants_analysis.R to visualize 




### generates reads for each base pair

# echo -e "chrom\pos\depth" > "$out_file"

# while read -r chr; do
#     echo "Processing $chr..."
#     samtools depth -r "$chr" "$bam" >> "$out_file"
# done < $chr_list

# tr '\t' ',' < /scratch/ejy4bu/compBio/bam_analysis/chr1_depth.tsv > /scratch/ejy4bu/compBio/bam_analysis/chr1_depth.csv
# rm -f $outfile


