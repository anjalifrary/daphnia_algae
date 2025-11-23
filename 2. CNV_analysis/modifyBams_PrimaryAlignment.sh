#!/usr/bin/env bash
#
#SBATCH -J cleanBams # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
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
bam=${bams[$SLURM_ARRAY_TASK_ID-1]}
samp=$(basename $bam .sort.bam)
samp_folder=${out_dir}/$samp
mkdir -p $samp_folder

### remove secondary 0x100 and supplementary 0x800 alignments
# new bam is copied to new folder in cnv 
samtools view -bh -F 0x900 $file -o ${samp_folder}/${samp}.clean.bam

### dedup cleaned bams 
echo "Deduplicating $samp"
dedup_bam=${samp_folder}/${samp}.dedup.bam
java -Xmx45G -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    I="$bam" \
    O="$dedup_bam" \
    M="${samp_folder}/${samp}.dedup.metrics" \
    REMOVE_DUPLICATES=true \
    CREATE_INDEX=true





# for file in $in_dir/*/*/*.sort.bam; do 
#     samp=$(basename $file .sort.bam)
#     echo "Copying $samp"
#     mkdir -p ${out_dir}/$samp
#     cp $file ${out_dir}/${samp}/${samp}.sort.bam

# done

# ### remove secondary 0x100 and supplementary 0x800 alignments
# for file in $out_dir/*/*/*.sort.bam; do
#     samp=$(basename $file .sort.bam)
#     echo "Processing $samp"
#     samp_folder=$(dirname "$file")
#     samtools view -bh -F 0x900 $file -o ${samp_folder}/${samp}.clean.bam
# done

# for bam in $out_dir/*/*/*.clean.bam; do
#     samp=$(basename $bam .clean.bam)
#     samp_folder=$(dirname "$bam")
#     dedup_bam=${samp_folder}/${samp}.dedup.bam
#     echo "Deduplicating $samp"
#     java -Xmx45G -jar $EBROOTPICARD/picard.jar MarkDuplicates \
#         I="$bam" \
#         O="$dedup_bam" \
#         M="${samp_folder}/${samp}.dedup.metrics" \
#         REMOVE_DUPLICATES=true \
#         CREATE_INDEX=true
# done



