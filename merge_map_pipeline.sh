#!/usr/bin/env bash
#SBATCH -J merge_pipeline    # Job name
#SBATCH --ntasks=1        # Single task per job
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH -N 1              # Run on one node
#SBATCH -t 0-10:00        # 10 hours runtime
#SBATCH --mem=100G        # Memory per node
#SBATCH -o /scratch/ejy4bu/erroroutputs/down.%A_%a.out  # Standard output
#SBATCH -e /scratch/ejy4bu/erroroutputs/down.%A_%a.err  # Standard error
#SBATCH -p standard       # Partition
#SBATCH --account=berglandlab

sbatch merge_fastq.sh
sbatch trim_fastq.sh
sbatch map_bam_ShortReads.sh
