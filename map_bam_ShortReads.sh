#!/usr/bin/env bash
#SBATCH -J makebams    # Job name
#SBATCH --ntasks=1        # Single task per job
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH -N 1              # Run on one node
#SBATCH -t 0-10:00        # 10 hours runtime
#SBATCH --mem=100G        # Memory per node
#SBATCH -o /scratch/ejy4bu/erroroutputs/map.%A_%a.out  # Standard output
#SBATCH -e /scratch/ejy4bu/erroroutputs/map.%A_%a.err  # Standard error
#SBATCH -p standard       # Partition
#SBATCH --account=berglandlab


#sbatch file.sh
#Submitted batch job 3436256
#sacct -j 3436256
#squeue -u (yourcomputingid)
#scancel 3436256 #to cancel the job

#ijob -A berglandlab -c10 -p standard --mem=40G



# Load necessary modules
module load gcc htslib
module load sratoolkit/3.1.1
#module load trimmomatic
module load bwa
module load samtools
#module load picard

:<<SRR
# Define working directories
infq="/scratch/ejy4bu/compBio/fastq"
outbam="/scratch/ejy4bu/compBio/bams"
#outbam="/scratch/ejy4bu/compBio/mapped_bam"

# Ensure output directories exist
mkdir -p "${infq}" "${outbam}" 
#"${outbam}"
SRR

# Extract fields (assuming CSV format: sample_id,reference_path)
ref_path=/project/berglandlab/chlorella_sequencing/reference_genome/GCA_023343905.1_cvul_genomic.fa


# array of sample directories for parallelization
#iterate through directories that contain the forward and reverse short read fastq files
#map to reference genome (assembled reads)
sample_folders=($(ls -d ${infq}/*/)) #Array to folder paths
cd "$infq"

samp_dir="${sample_folders[$SLURM_ARRAY_TASK_ID-1]}"
#samp="SRR14426882"
samp=$(basename "${samp_dir}")

echo "Entering folder: $samp"
cd "$samp" || continue  # Enter folder

#check for forward merged/trimmed fq.gz
if ls *trimmedmerged1.fq.gz 1> /dev/null 2>&1; then
    forward=(${samp}*trimmedmerged1.fq.gz)
else
    echo "Warning: No fastq in $samp"
fi
#check for reverse merged/trimmed fq.gz
if ls *trimmedmerged2.fq.gz 1> /dev/null 2>&1; then
    reverse=(${samp}*trimmedmerged2.fq.gz)
else
    echo "Warning: No fastq in $samp"
fi

echo "Processing sample : ${samp}"

bwa mem -t 10 -K 100000000 -Y ${ref_path} ${forward[@]} ${reverse[@]} | \
samtools view -uh -q 20 -F 0x100 | \
samtools sort --threads 10 -o "${outbam}/${samp}.sort.bam"

samtools index "${outbam}/${samp}.sort.bam"

:<<iterative
for samp_directory in ${infq}/*; 
    do
        samp=$(basename "${samp_directory}")

        forward=(${samp_directory}/*_1.fastq)
        reverse=(${samp_directory}/*_2.fastq)

        if [[ ! -f $forward || ! -f $reverse ]]; then
            echo "Skipping ${samp} (missing fastq files)"
            continue
        fi 

        echo "Processing sample : ${samp}"

        bwa mem -t 10 -K 100000000 -Y ${ref_path} ${forward} ${reverse} | \
        samtools view -uh -q 20 -F 0x100 | \
        samtools sort --threads 10 -o "${outbam}/${samp}.sort.bam"

        samtools index "${outbam}/${samp}.sort.bam"
done
iterative


:<<test
#test sample /scratch/ejy4bu/compBio/fastq/SRR14426881
# test="SRR14476638"
samp_directory="/scratch/ejy4bu/compBio/fastq/SRR14476638"
samp=$(basename "${samp_directory}")

forward=(${samp_directory}/*_1.fastq)
reverse=(${samp_directory}/*_2.fastq)

if [[ ! -f $forward || ! -f $reverse ]]; then
    echo "Skipping ${samp} (missing fastq files)"
fi 

echo "Processing sample : ${samp}"

bwa mem -t 10 -K 100000000 -Y ${ref_path} ${forward} ${reverse} | \
samtools view -uh -q 20 -F 0x100 | \
samtools sort --threads 10 -o "${outbam}/${samp}.sort.bam"

samtools index "${outbam}/${samp}.sort.bam"
test

:<<delete

forward=${samp_dir}/*_1.fastq
reverse=${samp_dir}/*_2.fastq

:<<skip_if
if [ ! -f $forward || ! -f $reverse ]; then
    echo "Skipping ${samp} (missing fastq files)"
fi 
skip_if
delete

# to make bam viewable as a sam file
# samtools view -h /scratch/ejy4bu/compBio/bams/SRR14426881.sort.bam > /scratch/ejy4bu/compBio/bams/SRR14426881.sort.sam
