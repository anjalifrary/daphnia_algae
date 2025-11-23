#!/usr/bin/env bash
#
#SBATCH -J cleanBams 
#SBATCH --array=1-380%15
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/err_outs/cleanBams/cleanBams.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/cleanBams/cleanBams.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load samtools
module load picard

### copy bams to new folder
in_dir="/scratch/ejy4bu/compBio/bams"
out_dir="/scratch/ejy4bu/compBio/cnv/bams"
mkdir -p $out_dir

bams=($in_dir/*/*/*.sort.bam)
bam=${bams[$((SLURM_ARRAY_TASK_ID-1))]}
samp=$(basename $bam .sort.bam)
samp_folder=${out_dir}/$samp
mkdir -p $samp_folder

### remove secondary 0x100 and supplementary 0x800 alignments
# new bam is copied to new folder in cnv 
clean_bam="${samp_folder}/${samp}.clean.bam"
if [ ! -f "$clean_bam" ]; then
    samtools view -bh -F 0x900 $bam -o $clean_bam
else
    echo "$clean_bam exists, skipping"
fi

# check number of reads in outfile
echo "Original reads: $(samtools view -c "$bam")"
echo "After cleaning: $(samtools view -c "$clean_bam")"

### dedup cleaned bams 
echo "Deduplicating $samp"
dedup_bam=${samp_folder}/${samp}.dedup.bam

# if [ -f "$dedup_bam" ]; then
#     echo "$dedup_bam already exists, skipping"
#     exit 0
# fi

java -Xmx45G -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    -I ${samp_folder}/${samp}.clean.bam \
    -O "$dedup_bam" \
    -M "${samp_folder}/${samp}.dedup.metrics" \
    -REMOVE_DUPLICATES true \
    -CREATE_INDEX true
