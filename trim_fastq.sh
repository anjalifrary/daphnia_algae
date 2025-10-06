#!/usr/bin/env bash
#
#SBATCH -J trimSRA # A single job name for the array
#SBATCH --ntasks-per-node=10 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/erroroutputs/trim.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/erroroutputs/trim.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

#ijob -A berglandlab -c10 -p standard --mem=40G

# Load necessary modules
module load gcc htslib
module load trimmomatic

#for pipeline, get sample_dir from array
sample_dir="$1"
samp_name=$(basename "$sample_dir")

for forward in "$sample_dir"/*_1.fq.gz; do
    reverse="${forward/_1.fq.gz/_2.fq.gz}"
    lane_name="${forward%_1.fq.gz}"

    echo "Processing lane: ${lane_name}"

    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 10 \
        ${lane_name}_1.fq.gz \
        ${lane_name}_2.fq.gz \
        ${lane_name}_1.P.trimm.fastq \
        ${lane_name}_1.U.trimm.fastq \
        ${lane_name}_2.P.trimm.fastq \
        ${lane_name}_2.U.trimm.fastq \
        ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10:8:true

    /home/ejy4bu/miniconda3/bin/pear \
        -f ${lane_name}_1.P.trimm.fastq \
        -r ${lane_name}_2.P.trimm.fastq \
        -o ${lane_name} \
        -j 10
    echo "Done working on $lane_name"

done
echo "All lanes trimmed for $samp_name"


