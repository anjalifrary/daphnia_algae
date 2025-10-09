#!/usr/bin/env bash
#SBATCH -J merge_pipeline    # Job name
#SBATCH --ntasks=1        # Single task per job
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH -N 1              # Run on one node
#SBATCH -t 0-10:00        # 10 hours runtime
#SBATCH --mem=100G        # Memory per node
#SBATCH -o /scratch/ejy4bu/erroroutputs/pipe.%A_%a.out  # Standard output
#SBATCH -e /scratch/ejy4bu/erroroutputs/pipe.%A_%a.err  # Standard error
#SBATCH -p standard       # Partition
#SBATCH --account=berglandlab

#ijob -A berglandlab -c10 -p standard --mem=100G


MY_DATA="/scratch/ejy4bu/compBio/fastq/Old_Algae_fastq"

for SAMPLE_DIR in "$MY_DATA"/*; do
    folder=$(basename "$SAMPLE_DIR")
    echo "Submitting jobs for sample: $folder"

    TRIM_JOB=$(sbatch --parsable trim_fastq.sh "$SAMPLE_DIR")

    MERGE_JOB=$(sbatch --parsable --dependency=afterok:$TRIM_JOB merge_fastq.sh "$SAMPLE_DIR")
done

MAP_JOB=$(sbatch --dependency=afterok:$MERGE_JOB map_bam_ShortReads.sh "$MY_DATA")


#SAMPLE="/scratch/ejy4bu/compBio/Robert_samples/RobertUK_F1"
#TRIM_JOB=$(sbatch --parsable trim_fastq.sh "$SAMPLE_DIR" "$SAMPLE")


