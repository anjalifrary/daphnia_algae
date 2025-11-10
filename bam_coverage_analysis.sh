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

module load samtools
module load bedtools

### 1. get coverage counts

#https://www.htslib.org/doc/samtools-depth.html

# samtools depth -flags in_bam.bam 

# https://www.htslib.org/doc/samtools-coverage.html

chr_list="/scratch/ejy4bu/compBio/genomefiles/ChrScaffoldList.txt"
bam="/scratch/ejy4bu/compBio/bams/Robert_samples_bams/RobertUK_F1/RobertUK_F1.dedup.bam"
out_file="/scratch/ejy4bu/compBio/bam_analysis/chr1_depth.tsv"

echo -e "chrom\pos\depth" > "$out_file"
samtools depth -r SIDB01000001.1 $bam >> $out_file



while read -r chr; do
    echo "Processing $chr..."
    samtools depth -r "$chr" "$bam" >> "$out_file"
done < $chr_list

tr '\t' ',' < /scratch/ejy4bu/compBio/bam_analysis/chr1_depth.tsv > /scratch/ejy4bu/compBio/bam_analysis/chr1_depth.csv
rm -f $outfile


### 2. normalize coverage, using R

### 3. 
