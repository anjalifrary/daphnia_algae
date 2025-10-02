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

:<<copy_files
DATA_DIRS=(
"/project/berglandlab/Robert/shortread_data/data/01.RawData/RobertUK_F1"
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
mkdir -p "/scratch/ejy4bu/compBio/Robert_samples"
for samp in "${!DATA_DIRS[@]}"; do
    samp_dir=${DATA_DIRS[$samp]}
    name=$(basename "$samp_dir")
    mkdir -p "/scratch/ejy4bu/compBio/Robert_samples/$name"
    files=("$samp_dir"/*.fq.gz)
    if [ ${#files[@]} -eq 0 ]; then
        echo "No .fq.gz files in $samp_dir, skipping"
        continue
    fi
    cp -n "${files[@]}" "/scratch/ejy4bu/compBio/Robert_samples/${name}"
    echo "Copied ${#files[@]} files"
done
copy_files

MY_DATA="/scratch/ejy4bu/compBio/Robert_samples"
cd "$MY_DATA"

for SAMPLE_DIR in "$MY_DATA"/*; do
    SAMPLE=$(basename "$SAMPLE_DIR")
    echo "Submitting jobs for sample: $SAMPLE"
    cd "$SAMPLE_DIR" || continue
    gunzip *.fq.gz
    echo "Unzipping fq.gz files: $SAMPLE"
    #trim
    TRIM_JOB=$(sbatch --parsable trim_fastq.sh "$SAMPLE_DIR" "$SAMPLE")
    #merge
    MERGE_JOB=$(sbatch --parsable --dependency=afterok:$TRIM_JOB merge_fastq.sh "$SAMPLE_DIR" "$SAMPLE")
    #map
    MAP_JOB=$(sbatch --dependency=afterok:$MERGE_JOB map_bam_ShortReads.sh "$SAMPLE_DIR" "$SAMPLE")
done

#SAMPLE="/scratch/ejy4bu/compBio/Robert_samples/RobertUK_F1"
#TRIM_JOB=$(sbatch --parsable trim_fastq.sh "$SAMPLE_DIR" "$SAMPLE")


