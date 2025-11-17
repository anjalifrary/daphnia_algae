#!/usr/bin/env bash
#SBATCH -J convertVCF    # Job name
#SBATCH --ntasks=1        # Single task per job
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH -N 1              # Run on one node
#SBATCH -t 0-10:00        # 10 hours runtime
#SBATCH --mem=50G        # Memory per node
#SBATCH -o /scratch/ejy4bu/erroroutputs/convert.%A_%a.out  # Standard output
#SBATCH -e /scratch/ejy4bu/erroroutputs/convert.%A_%a.err  # Standard error
#SBATCH -p standard       # Partition
#SBATCH --account=berglandlab

module load gatk

vcf="/project/berglandlab/chlorella_sequencing/combined_chlorella_annotated.vcf"
outfile="/project/berglandlab/chlorella_sequencing/combined_chlorella_annotated.tsv"

gatk VariantsToTable \
-V $vcf \
-O $outfile \
-F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F INFO \
-GF GT 

#-GF AD -GF DP -GF GQ
