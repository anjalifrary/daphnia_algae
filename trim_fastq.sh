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

# Load necessary modules
module load gcc htslib
# module load sratoolkit/3.1.1
module load trimmomatic
# module load bwa
# module load samtools
#module load picard

:<<download_pear
download miniconda and then pear: run in SSH terminal: (it will take a long time)
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
    source $HOME/miniconda3/bin/activate
    conda install -c bioconda pear

download_pear

#mv /scratch/ejy4bu/UK2022_2024/allshortreads/01.RawData/SR* /scratch/rjp5nc/UK2022_2024/allshortreads/01.RawData/SRR/

infq="/scratch/ejy4bu/compBio/fastq"

#PARENT_DIR="/scratch/ejy4bu/UK2022_2024/allshortreads/01.RawData/SRR"


# Get list of all folders
sample_folders=($(ls -d ${infq}/*/))  # Array of folder paths
# Select the folder based on the job array index
sample_dir="${sample_folders[$SLURM_ARRAY_TASK_ID - 1]}"
#sample_dir="/scratch/ejy4bu/compBio/fastq/SRR14476638"
# Extract folder name
samp_name=$(basename "$sample_dir")

# Change to the working directory
cd "$sample_dir" || exit 1


# Check if fastq files exist
#https://www.baeldung.com/linux/compgen-command-usage
if compgen -f -G "*.fastq" > /dev/null 2>&1; then
    echo "Processing $samp_name"
    # Run Trimmomatic
    trimmomatic PE -threads 10 \
        ${samp_name}_1.fastq \
        ${samp_name}_2.fastq \
        ${samp_name}_1.P.trimm.fastq \
        ${samp_name}_1.U.trimm.fastq \
        ${samp_name}_2.P.trimm.fastq \
        ${samp_name}_2.U.trimm.fastq \
        ILLUMINACLIP:/home/ejy4bu/miniconda3/bin/trimmomatic/TrimmomaticAdaptors/CombinedPE-PE.fa:2:30:10:8:true
    # Run PEAR to merge overlapping reads
    /home/ejy4bu/miniconda3/bin/pear \
        -f ${samp_name}_1.P.trimm.fastq \
        -r ${samp_name}_2.P.trimm.fastq \
        -o ${samp_name} \
        -j 10
else
    echo "Warning: No fastq files in $samp_name"
fi

:<<deletes

forward=${samp_dir}/*_1.fastq
reverse=${samp_dir}/*_2.fastq

if ![ -f $forward || -f $reverse ]; then
    echo "Skipping ${samp} (missing fastq files)"
fi 


    trimmomatic PE -threads 10 \ #replaced by ##deleted because installing trimmomatic
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 10 \
    ILLUMINACLIP:/home/ejy4bu/miniconda3/bin/trimmomatic/TrimmomaticAdaptors/CombinedPE-PE.fa:2:30:10:8:true
        
    #deleted because downloading trimmomatic instead
        ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10:8:true

    /home/ejy4bu/miniconda3/bin/pear \

    # could do:

    export EBROOT_PEAR=/home/ejy4bu/miniconda3
    $EBROOT_PEAR/bin/pear -f sample_1.P.trimm.fastq -r sample_2.P.trimm.fastq -o sample -j 10

    #didn't work:
    $CONDA_PREFIX/bin/pear -f sample_1.P.trimm.fastq -r sample_2.P.trimm.fastq -o sample -j 10

    #no pear module on the cluster
        pear -f sample_1.P.trimm.fastq -r sample_2.P.trimm.fastq -o sample -j 10


    #calling PEAR module. 
        -f ${samp_name}_1.P.trimm.fastq \
        -r ${samp_name}_2.P.trimm.fastq \
        -o ${samp_name} \
        -j 10



deletes