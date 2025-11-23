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
module load blast


proj="/project/berglandlab/chlorella_sequencing/reference_genome/GCA_023343905.1_cvul_genomic.fa"
out_dir="/scratch/ejy4bu/compBio/cnv/reference"
mkdir -p $out_dir

# cp $proj $out_dir
ref="${out_dir}/GCA_023343905.1_cvul_genomic.fa"
cleaned="${out_dir}/GCA_023343905.1_cvul_genomic.cleaned.fasta"

# cd $out_dir

seqtk seq $ref > $cleaned

#make blast db from cleaned reference genome (type nucleic acid)
makeblastdb -in $cleaned -out chlorellaDB -dbtype nucl -title chlorellaDB -parse_seqids

RepeatModeler -database chlorellaDB -pa 10 -engine ncbi