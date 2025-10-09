#!/usr/bin/env bash
#SBATCH -J map_pipeline    # Job name
#SBATCH --ntasks=1        # Single task per job
#SBATCH --cpus-per-task=4 # Number of CPU cores per task
#SBATCH -N 1              # Run on one node
#SBATCH -t 0-10:00        # 10 hours runtime
#SBATCH --mem=50G        # Memory per node
#SBATCH -o /scratch/ejy4bu/erroroutputs/pipe.%A_%a.out  # Standard output
#SBATCH -e /scratch/ejy4bu/erroroutputs/pipe.%A_%a.err  # Standard error
#SBATCH -p standard       # Partition
#SBATCH --account=berglandlab

#after running the trim and merge pipeline successfully: 

MY_DATA="/scratch/ejy4bu/compBio/Robert_samples"

MAP_JOB=$(sbatch --parsable map_bam_ShortReads.sh "$MY_DATA")
