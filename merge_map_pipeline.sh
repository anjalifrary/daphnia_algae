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

DATA_DIRS=(
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_F1"
)
:<<theRest
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_F2"
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_F3"
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_F4"
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_F5"
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_F6"
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_F7"
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_F8"
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_F9"
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_F10"
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_F11"
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_F12"
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_G1"
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_G2"
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_G3"
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_G4"
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_G5"
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_G6"
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_G7"
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_G8"
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_G9"
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_G10"
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_G11"
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_G12"
)
theRest

mkdir -p "/scratch/ejy4bu/compBio/Robert_samples"
#my_dir="/scratch/ejy4bu/compBio/Robert_samples"

for samp in "${!DATA_DIRS[@]}"; do
    samp_dir=${DATA_DIRS[$samp]}
    name=$(basename "$samp_dir")
    mkdir -p "/scratch/ejy4bu/compBio/Robert_samples/$samp"
    files=("$samp_dir"/*.fq.gz)

    cp -n "${files[@]}/" "/scratch/ejy4bu/compBio/Robert_samples/${name}"
done


:<<later
for SAMPLE_DIR in "${DATA_DIRS[@]}"; do
    SAMPLE=$(basename "$SAMPLE_DIR")
    
    #trim
    TRIM_JOB=$(sbatch --parsable --dependency=afterok:$MERGE_JOB trim_fastq.sh "$SAMPLE" )

    #merge
    MERGE_JOB=$(sbatch --parsable merge_fastq.sh "$SAMPLE_DIR" "$SAMPLE")

    #map
    sbatch --dependency=afterok:$TRIM_JOB map_bam_ShortReads.sh "$SAMPLE"
done


sbatch merge_fastq.sh
sbatch trim_fastq.sh
sbatch map_bam_ShortReads.sh
later
