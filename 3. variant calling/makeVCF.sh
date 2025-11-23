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
ref_fasta=/project/berglandlab/chlorella_sequencing/reference_genome/GCA_023343905.1_cvul_genomic.fa
out_vcf=/scratch/ejy4bu/compBio/vcfs
chr_list=/scratch/ejy4bu/compBio/genomefiles/ChrScaffoldList.txt

mkdir -p $out_vcf

# Get the chromosome for the current SLURM task


SLURM_ARRAY_TASK_ID=1

chr="SIDB01000002.1"
#chr=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $chr_list)
echo $chr

bam_list="/scratch/ejy4bu/compBio/vcfs/bam_list.txt"
ls $bam_root/*/*/*.dedup.bam > $bam_list

# # Store BAM file list in a variable to maintain order
# bam_list=$(ls $bam_root/*/*/*.dedup.bam | sort)
# echo "Using ${#bam_list[@]} BAMs"

# make sure sample names are specified
    # flag in bam file stores sample name, make sure it is passed to vcf
    # see how mpileup or varscan work

# Generate VCF with VarScan
samtools mpileup \
    -r ${chr} \
    --fasta-ref $ref_fasta \
    -b $bam_list| \
java -jar $EBROOTVARSCAN/VarScan.v2.4.4.jar mpileup2snp \
    /dev/stdin \
    --min-coverage 4 \
    --min-var-freq 0.001 \
    --output-vcf > $out_vcf/${chr}.vcf

# for vcf in /project/berglandlab/chlorella_sequencing/vcfs/*.vcf; do
#     cp $vcf /scratch/ejy4bu/compBio/vcfs;
# done