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


MY_DATA="/scratch/ejy4bu/compBio/Robert_samples"
#cd "$MY_DATA"

#SAMPLE_DIR="/scratch/ejy4bu/compBio/Robert_samples/RobertUK_F1"

for SAMPLE_DIR in "$MY_DATA"/*; do
    SAMPLE=$(basename "$SAMPLE_DIR")
    echo "Submitting jobs for sample: $SAMPLE"
    #cd "$SAMPLE_DIR" || continue
    #gunzip *.fq.gz
    #echo "Unzipping fq.gz files: $SAMPLE"
    #trim
    for samp in "$SAMPLE_DIR"; do
    TRIM_JOB=$(sbatch --parsable trim_fastq.sh "$SAMPLE_DIR" "$SAMPLE")
    done
    #merge
    MERGE_JOB=$(sbatch --parsable --dependency=afterok:$TRIM_JOB merge_fastq.sh "$SAMPLE_DIR" "$SAMPLE")
    #map
    MAP_JOB=$(sbatch --dependency=afterok:$MERGE_JOB map_bam_ShortReads.sh "$SAMPLE_DIR" "$SAMPLE")
done

#SAMPLE="/scratch/ejy4bu/compBio/Robert_samples/RobertUK_F1"
#TRIM_JOB=$(sbatch --parsable trim_fastq.sh "$SAMPLE_DIR" "$SAMPLE")


