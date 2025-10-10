#!/usr/bin/env bash
#SBATCH -J makebams    # Job name
#SBATCH --array=1-320
#SBATCH --ntasks=1        # Single task per job
#SBATCH --cpus-per-task=4 # Number of CPU cores per task
#SBATCH -N 1              # Run on one node
#SBATCH -t 0-10:00        # 10 hours runtime
#SBATCH --mem=50G        # Memory per node
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
module load bwa
module load samtools

# Extract fields (assuming CSV format: sample_id,reference_path)
ref_path=/project/berglandlab/chlorella_sequencing/reference_genome/GCA_023343905.1_cvul_genomic.fa
#infq="$1"
infq="/scratch/ejy4bu/compBio/fastq/Old_Algae_fastq"

# array of sample directories for parallelization
#iterate through directories that contain the forward and reverse short read fastq files
#map to reference genome (assembled reads)

sample_folders=($(ls -d ${infq}/*/)) #Array to folder paths
samp_dir="${sample_folders[$SLURM_ARRAY_TASK_ID-1]}"
#samp_dir="/scratch/ejy4bu/compBio/Robert_samples/RobertUK_F1"
samp=$(basename "${samp_dir}")

outbam="/scratch/ejy4bu/compBio/bams/Old_Algae_bams/${samp}"
mkdir -p "${outbam}" 

if ls ${outbam}/*.sort.bam 1> /dev/null 2>&1; then
    echo "Already mapped ${samp}"
    exit 1
fi

#check for forward merged/trimmed fq.gz
if ls ${samp_dir}/*_trimmedmerged1.fq.gz 1> /dev/null 2>&1; then
    forward=$(ls ${samp_dir}/*_trimmedmerged1.fq.gz)
else
    echo "Warning: No fastq in $samp"
fi
#check for reverse merged/trimmed fq.gz
if ls ${samp_dir}/*trimmedmerged2.fq.gz 1> /dev/null 2>&1; then
    reverse=$(ls ${samp_dir}/*_trimmedmerged2.fq.gz)
else
    echo "Warning: No fastq in $samp"
fi

echo "Processing sample : ${samp}"

bwa mem -t 10 -K 100000000 -Y "$ref_path" "$forward" "$reverse" | \
samtools view -uh -q 20 | \
samtools sort --threads 10 -o "$outbam/${samp}.sort.bam"

samtools index "$outbam/${samp}.sort.bam"
echo "finished mapping $samp"


# to make bam viewable as a sam file
# samtools view -h /scratch/ejy4bu/compBio/bams/SRR14426881.sort.bam > /scratch/ejy4bu/compBio/bams/SRR14426881.sort.sam
# samtools idxstats
    # check for reads in * 0 (last line)
    # sum values for total reads
    # sum reads in top 14
#pos and neg controls for chlorella

#samtools make a pileup file 
#turn that into vcf file
# keep counts of alt depth and ref depth 