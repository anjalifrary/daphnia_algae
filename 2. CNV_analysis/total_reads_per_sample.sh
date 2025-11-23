#!/usr/bin/env bash
#
#SBATCH -J subset 
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH -N 1 # on one node
#SBATCH -t 0-10:00 # 10 hours
#SBATCH --mem 100G
#SBATCH -o /scratch/ejy4bu/err_outs/subset/subset.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/subset/subset.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

module load samtools

bam_dir="/scratch/ejy4bu/compBio/cnv/bams"
out_file="/scratch/ejy4bu/compBio/cnv/bam_read_counts.csv"

# write header to csv output file
echo "sample,group,total_reads" > $out_file

for bam in ${bam_dir}/*/*/*.dedup.bam; do
    sample=$(basename $bam .dedup.bam)
    group=$(basename $(dirname $(dirname $bam)))

    echo "Processing sample $sample in group $group"

    # count total mapped reads, -F 4 exclude unmapped reads
    reads=$(samtools view -c -F 4 "$bam")
    echo "${sample},${group},${reads}" >> $out_file
done

echo "Finished calculating reads for all samples"


### getting text list of sample subsets

# infile="/scratch/ejy4bu/compBio/cnv/subset_samples.csv"
# outfile="/scratch/ejy4bu/compBio/cnv/subset_sample_list.txt"

## include all samples
# tail -n +2 $infile | cut -d',' -f1 > $outfile

## remove control group
# tail -n +2 $infile | awk -F, '$3 != "REED_Sephedex" {print $1}' > $outfile

### path names for each subset of samples

:<<path_names
bam_dir="/scratch/ejy4bu/compBio/cnv/bams"
sample_list="/scratch/ejy4bu/compBio/cnv/data_tables/subset_samples.csv"
targ_group="UTEX"
bam_list="/scratch/ejy4bu/compBio/cnv/data_tables/${targ_group}_bamList.txt"

> $bam_list

tail -n +2 $sample_list | while IFS=',' read -r sample Date group total_reads sum_total_reads; do
    sample=$(echo "$sample" | tr -d '[] \r')
    group=$(echo "$group" | tr -d '[] \r')
    if [[ "$group" == "$targ_group" ]]; then
        echo ${bam_dir}/${group}/${sample}/${sample}.dedup.bam >> $bam_list
    fi
done
path_names




