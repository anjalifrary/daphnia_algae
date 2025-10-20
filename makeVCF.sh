#!/usr/bin/env bash
#
#SBATCH -J makeVCF # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/erroroutputs/makeVCF/vcf.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/erroroutputs/makeVCF/vcf.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

mkdir -p /scratch/ejy4bu/erroroutputs/makeVCF

module load samtools varscan bcftools

# sbatch --array=1-14 ~/daphnia_algae/makeVCF.sh

bam_root=/scratch/ejy4bu/compBio/bams
ref_fasta=/project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.fa
out_vcf=/scratch/ejy4bu/compBio/vcfs

mkdir -p $out_vcf

# Get the chromosome for the current SLURM task
chr=$(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/ejy4bu/compBio/genomefiles/ChrScaffoldList)
echo $chr

# Store BAM file list in a variable to maintain order
bam_list=$(ls $bam_root/*/*/*.dedup.bam | sort)
echo "Using ${#bam_list[@]} BAMs"

# Generate VCF with VarScan
samtools mpileup -@ 10 \
    -r ${chr} \
    --fasta-ref $ref_fasta \
    $bam_list | \
java -jar $EBROOTVARSCAN/VarScan.v2.4.4.jar mpileup2snp \
    /dev/stdin \
    --min-coverage 1 \
    --min-var-freq 0.001 \
    --output-vcf > $out_vcf/${chr}.vcf


# I lowered the min-coverage to 1 (originally 4) because my vcf files were coming up empty...
