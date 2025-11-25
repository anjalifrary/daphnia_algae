#!/usr/bin/env bash
#
#SBATCH -J repeatModeler 
#SBATCH --cpus-per-task=10
#SBATCH --ntasks=1
#SBATCH -N 1 # on one node
#SBATCH -t 5-00:00 # 5 days
#SBATCH --mem 150G
#SBATCH -o /scratch/ejy4bu/err_outs/repeatModeler/repeatModeler.%A_%a.out # Standard output
#SBATCH -e /scratch/ejy4bu/err_outs/repeatModeler/repeatModeler.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

# make conda environment with repeatmodeler in it

#conda create -n repeatmodeler_env -c bioconda -c conda-forge repeatmodeler perl=5.22.0
source ~/miniconda3/etc/profile.d/conda.sh

conda activate repeatmodeler_env

# proj="/project/berglandlab/chlorella_sequencing/reference_genome/GCA_023343905.1_cvul_genomic.fa"
out_dir="/scratch/ejy4bu/compBio/cnv/reference"
# mkdir -p $out_dir

# cp $proj $out_dir
# ref="${out_dir}/GCA_023343905.1_cvul_genomic.fa"
cleaned="${out_dir}/GCA_023343905.1_cvul_genomic.cleaned.fasta"
# seqtk seq $ref > $cleaned

cd $out_dir

ls -lh $cleaned


#make blast db from cleaned reference genome (type nucleic acid)
makeblastdb -in $cleaned -out my_db -dbtype nucl -title my_db -parse_seqids
BuildDatabase -name my_db $cleaned

RepeatModeler -database my_db -pa 10



# make index and produce scaffold length file
out_dir="/scratch/ejy4bu/compBio/cnv/reference"
cleaned="${out_dir}/GCA_023343905.1_cvul_genomic.cleaned.fasta"

module load samtools
mkdir -p "scratch/ejy4bu/compBio/cnv/megabams/debug"
scaffold_lengths="scratch/ejy4bu/compBio/cnv/megabams/debug/scaffold_lengths.txt"
samtools faidx "$cleaned"

# use the .fai to get scaffold names + lengths
fai="${cleaned}.fai"

cut -f1,2 "$fai" > "$scaffold_lengths"

# quick checks
wc -l "$scaffold_lengths"            # number of scaffolds
head -n 10 "$scaffold_lengths"       # first 10 scaffolds with lengths
tail -n 10 "$scaffold_lengths"       # last 10 scaffolds with lengths

# total bases
awk '{s+=$2} END{print "total_bp="s}' "$scaffold_lengths"

cut -f2 $fai | sort -nr > "scratch/ejy4bu/compBio/cnv/megabams/debug/lengths.txt"
awk 'BEGIN{sum=0; i=0} {sum+=$1; a[i++]=$1; total+= $1} END{half=total/2; s=0; for(j=0;j<i;j++){s+=a[j]; if(s>=half){print "N50="a[j]; break}} }' "scratch/ejy4bu/compBio/cnv/megabams/debug/lengths.txt"

consensi="/sfs/weka/scratch/ejy4bu/compBio/cnv/reference/RM_175000.SunNov231803172025/consensi.fa"

# change path to the actual consensi file found above
grep -c '^>' $consensi
grep -c '^>' $consensi.classified   # if present
