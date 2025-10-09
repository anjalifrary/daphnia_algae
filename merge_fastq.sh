#!/usr/bin/env bash
#
#SBATCH -J merge # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/erroroutputs/merge.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/erroroutputs/merge.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load necessary modules

# Set the parent directory containing all the folders
# PARENT_DIR="/scratch/ejy4bu/UK2022_2024/allshortreads/01.RawData/newseq/"  # Change this to your actual path

#get sample from pipeline 
sample_dir="$1"
samp_name=$(basename"${sample_dir}")

if ls "$sample_dir"/*_1.P.trimm.fastq 1> /dev/null 2>&1; then
    echo "Merging all forward reads in $samp_name..."
    cat "$sample_dir"/*_1.P.trimm.fastq | gzip > "$sample_dir/${samp_name}_trimmedmerged1.fq.gz"
    rm "$sample_dir"/*_1.P.trimm.fastq
    rm "$sample_dir"/*_1.fastq
else
    echo "Warning: No fastq in $samp_name"
fi

if ls "$sample_dir"/*_2.P.trimm.fastq 1> /dev/null 2>&1; then
    echo "Merging all reverse reads in $samp_name..."
    cat "$sample_dir"/*_2.P.trimm.fastq | gzip > "$sample_dir/${samp_name}_trimmedmerged2.fq.gz"
    rm "$sample_dir"/*_2.P.trimm.fastq
    rm "$sample_dir"/*_2.fastq
else
    echo "Warning: No reverse reads found in $samp_name"
fi

echo "done merging in $samp_name"