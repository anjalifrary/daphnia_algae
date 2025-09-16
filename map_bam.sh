#!/usr/bin/env bash
#SBATCH -J makebams    # Job name
#SBATCH --ntasks=1        # Single task per job
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH -N 1              # Run on one node
#SBATCH -t 0-10:00        # 10 hours runtime
#SBATCH --mem=100G        # Memory per node
#SBATCH -o /scratch/ejy4bu/erroroutputs/down.%A_%a.out  # Standard output
#SBATCH -e /scratch/ejy4bu/erroroutputs/down.%A_%a.err  # Standard error
#SBATCH -p standard       # Partition
#SBATCH --account=berglandlab


#sbatch file.sh
#Submitted batch job 3436256
#sacct -j 3436256
#squeue -u (yourcomputingid)
#scancel 3436256 #to cancel the job





# Load necessary modules
module load gcc htslib
module load sratoolkit/3.1.1
module load trimmomatic
module load bwa
module load samtools
module load picard

# Define working directories
infq="/scratch/ejy4bu/compBio/fastq"
outfq="/scratch/ejy4bu/compBio/bams"
#outbam="/scratch/ejy4bu/compBio/mapped_bam"

# Ensure output directories exist
mkdir -p "${infq}" "${outfq}" "${outbam}"

# Extract fields (assuming CSV format: sample_id,reference_path)
ref_path=/project/berglandlab/chlorella_sequencing/reference_genome/GCA_023343905.1_cvul_genomic.fa


#ijob -A berglandlab -c10 -p standard --mem=40G

samp=long_read_Chlorella_read

# Map to reference genome (assembled reads)
bwa mem -t 10 -K 100000000 -Y ${ref_path} /project/berglandlab/chlorella_sequencing/raw_longread_from_Reed/m84128_250121_222443_s2.hifi_reads.bc2104.fq.gz |
samtools view -uh -q 20 -F 0x100 | \
samtools sort --threads 10 -o ${outfq}/${samp}.sort.bam

samtools index ${outfq}/${samp}.sort.bam

samtools view -h /scratch/ejy4bu/compBio/bams/chlorella_Reed.sort.bam > /scratch/ejy4bu/compBio/bams/chlorella_Reed.sort.sam

# /project/berglandlab/chlorella_sequencing/HMW/HMWDNAElvis3/m84128_250121_222443_s2.hifi_reads.bc2104.fastq | \
#-F 0x100 is to map secondary reads (repetitive regions)
