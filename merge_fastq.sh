#!/usr/bin/env bash
#
#SBATCH -J merge # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 2-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/erroroutputs/merge.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/erroroutputs/merge.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# Load necessary modules

# Set the parent directory containing all the folders
# PARENT_DIR="/scratch/ejy4bu/UK2022_2024/allshortreads/01.RawData/newseq/"  # Change this to your actual path

#get sample from array 
sample_dir="$1"
samp_name=$(basename "$sample_dir")
cd "$sample_dir" || exit 1



:<<SRR_samples
infq="/scratch/ejy4bu/compBio/fastq"
#sample_folders="/scratch/ejy4bu/compBio/fastq"
cd "$sample_folders" || exit
# Loop through all subdirectories
for folder in ${infq}/*; do
    samp=$(basename "${folder%/}")
    # samp="SRR14426882"
    echo "Entering folder: $samp"
    cd "$samp" || continue  # Enter folder
    # Check for forward reads
    if ls *_1.P.trimm.fastq 1> /dev/null 2>&1; then
        echo "Merging all forward reads in $samp..."
        cat *_1.P.trimm.fastq | gzip > "${samp}_trimmedmerged1.fq.gz"
        rm *_1.P.trimm.fastq
        rm *_1.fastq
    else
        echo "Warning: No fastq in $samp"
    fi
    # Check for reverse reads
    if ls *_2.P.trimm.fastq 1> /dev/null 2>&1; then
        echo "Merging all reverse reads in $samp..."
        cat *_2.P.trimm.fastq | gzip > "${samp}_trimmedmerged2.fq.gz"
        rm *_2.P.trimm.fastq
        rm *_2.fastq
    else
        echo "Warning: No reverse reads found in $samp"
    fi
    cd ..  # Move back to parent directory
done
echo "Done merging"
SRR_samples
