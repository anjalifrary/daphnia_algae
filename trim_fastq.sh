#!/usr/bin/env bash
#
#SBATCH -J trimSRA # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/erroroutputs/down.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/erroroutputs/down.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab_standard
#SBATCH --array=1-21%10   # Adjust the range based on the number of folders


#mv /scratch/ejy4bu/UK2022_2024/allshortreads/01.RawData/SR* /scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/SRR/

infq="/scratch/ejy4bu/compBio/fastq"

#PARENT_DIR="/scratch/ejy4bu/UK2022_2024/allshortreads/01.RawData/SRR"


# Get list of all folders
sample_folders=($(ls -d ${infq}/*/))  # Array of folder paths
# Select the folder based on the job array index
sample_dir="${sample_folders[$SLURM_ARRAY_TASK_ID - 1]}"
# Extract folder name
samp_name=$(basename "$sample_dir")

# Change to the working directory
cd "$sample_dir" || exit 1

# Check if fastq.gz files exist
if ls *.fastq.gz 1> /dev/null 2>&1; then
    echo "Processing $FOLDER_NAME"
    # Run Trimmomatic
    trimmomatic PE -threads 10 \
        ${FOLDER_NAME}_1.fastq.gz \
        ${FOLDER_NAME}_2.fastq.gz \
        ${FOLDER_NAME}_1.P.trimm.fastq.gz \
        ${FOLDER_NAME}_1.U.trimm.fastq.gz \
        ${FOLDER_NAME}_2.P.trimm.fastq.gz \
        ${FOLDER_NAME}_2.U.trimm.fastq.gz \
        ILLUMINACLIP:/home/rjp5nc/miniconda3/bin/trimmomatic/TrimmomaticAdaptors/CombinedPE-PE.fa:2:30:10:8:true
    # Run PEAR to merge overlapping reads
    /home/rjp5nc/miniconda3/bin/pear \
        -f ${FOLDER_NAME}_1.P.trimm.fastq.gz \
        -r ${FOLDER_NAME}_2.P.trimm.fastq.gz \
        -o ${FOLDER_NAME} \
        -j 10
else
    echo "Warning: No fastq.gz files in $FOLDER_NAME"
fi