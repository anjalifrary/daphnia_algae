#!/usr/bin/env bash
#SBATCH -J merge_pipeline_array    # Job name
#SBATCH --array=0-1
#SBATCH --ntasks=1        # Single task per job
#SBATCH --cpus-per-task=10 # Number of CPU cores per task
#SBATCH -N 1              # Run on one node
#SBATCH -t 0-10:00        # 10 hours runtime
#SBATCH --mem=50G        # Memory per node
#SBATCH -o /scratch/ejy4bu/erroroutputs/pipe.%A_%a.out  # Standard output
#SBATCH -e /scratch/ejy4bu/erroroutputs/pipe.%A_%a.err  # Standard error
#SBATCH -p standard       # Partition
#SBATCH --account=berglandlab
#SBATCH --begin=now+1hour

#ijob -A berglandlab -c10 -p standard --mem=100G


MY_DATA="/scratch/ejy4bu/compBio/fastq/Old_Algae_fastq"
# TRIM_SCRIPT="/home/ejy4bu/daphnia_algae/trim_fastq.sh"
# MERGE_SCRIPT="/home/ejy4bu/daphnia_algae/merge_fastq.sh"
# MAP_SCRIPT="/home/ejy4bu/daphnia_algae/map_bam_ShortReads.sh"

# chmod +x "$TRIM_SCRIPT"
# chmod +x "$MERGE_SCRIPT"
# chmod +x "$MAP_SCRIPT"

SAMPLES=($(ls -d ${MY_DATA}/*/))

echo "Samples = ${SAMPLES}"

SAMPLE_DIR="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"
samp_name=$(basename "$SAMPLE_DIR")

echo "Processing sample ${samp_name}. (Array task ID: $SLURM_ARRAY_TASK_ID)"

TRIM_JOB=$(sbatch --parsable trim_fastq.sh "$SAMPLE_DIR")

MERGE_JOB=$(sbatch --parsable --dependency=afterok:$TRIM_JOB merge_fastq.sh "$SAMPLE_DIR")

MAP_JOB=$(sbatch --parsable --dependency=afterok:$MERGE_JOB map_bam_ShortReads.sh "$MY_DATA")


echo "Submitted jobs for $samp_name:"
echo "  TRIM_JOB = $TRIM_JOB"
echo "  MERGE_JOB = $MERGE_JOB"
echo "  MAP_JOB = $MAP_JOB"
