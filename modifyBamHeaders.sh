#!/usr/bin/env bash
#SBATCH -J addRG_inplace
#SBATCH --ntasks=1
#SBATCH -N 1
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem=50G
#SBATCH -o /scratch/ejy4bu/erroroutputs/modifyBams/addRG_inplace.%A.out
#SBATCH -e /scratch/ejy4bu/erroroutputs/modifyBams/addRG_inplace.%A.err
#SBATCH -p standard
#SBATCH --account berglandlab

mkdir -p /scratch/ejy4bu/erroroutputs/modifyBams

module load samtools

bam_root=/scratch/ejy4bu/compBio/bams

for f in $(find "$bam_root" -type f -name "*.dedup.bam"); do
    sample=$(basename "$f" .dedup.bam)
    echo "Adding read group to $sample"
    
    tmp="${f%.bam}.tmp.bam"
    
    samtools addreplacerg \
        -r "ID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
        -o "$tmp" "$f"

    if [ $? -eq 0 ]; then
        mv "$tmp" "$f"
        samtools index "$f"
    else
        echo "Failed on $f"
        rm -f "$tmp"
    fi
done

# # extract chromosome lengths
# genome_lengths="/scratch/ejy4bu/compBio/genome_files/scaffold_lengths.txt"
# samtools faidx /project/berglandlab/chlorella_sequencing/reference_genome/GCA_023343905.1_cvul_genomic.fa
# cut -f1,2 /project/berglandlab/chlorella_sequencing/reference_genome/GCA_023343905.1_cvul_genomic.fa.fai > "$genome_lengths"
