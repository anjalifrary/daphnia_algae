#!/usr/bin/env bash
#
#SBATCH -J download_SRA       # A single job name for the array
#SBATCH --ntasks-per-node=10  # one core
#SBATCH -N 1                  # on one node
#SBATCH -t 6:00:00 ### 6 hours
#SBATCH --mem 10G
#SBATCH -o /scratch/ejy4bu/compBio/logs/prefetch.%A_%a.out # Standard output %A_%a contains job info 
#SBATCH -e /scratch/ejy4bu/compBio/logs/prefetch.%A_%a.err # Standard error
#SBATCH -p standard
#SBATCH --account berglandlab

#sacct -j 3283279 // to see the updated job status


wd=/scratch/ejy4bu/compBio
### run as: sbatch --array=1-$( wc -l < /project/berglandlab/anjali/metadata/algae_paths_anjali.csv )%10 ~/daphnia_algae/downloadSRA.sh
### sacct -j 64052181
### cat /scratch/ejy4bu/compBio/logs/prefetch.52222298_*.out | grep -B1 "do not"
### cat /scratch/ejy4bu/compBio/logs/prefetch.52222298_52.out

module load gcc/11.4.0 sratoolkit/3.1.1 

# SLURM_ARRAY_TASK_ID=1

CSV="/project/berglandlab/anjali/metadata/algae_paths_anjali.csv"

sranum=$(tail -n +2 "$CSV" | sed "${SLURM_ARRAY_TASK_ID}q;d" | cut -f8 -d",")
sampName=$sranum
# proj=$( sed "${SLURM_ARRAY_TASK_ID}q;d" "$CSV" | cut -f2 -d',' )

echo $sampName " / " $sranum 
# " / " $proj

### sranum=SRR1184609; proj=PRJNA194129

if [ ! -d "/scratch/ejy4bu/compBio/fastq/Old_Algae_fastq" ]; then
  mkdir /scratch/ejy4bu/compBio/fastq/Old_Algae_fastq
fi

if [ ! -d "/scratch/ejy4bu/compBio/fastq/Old_Algae_fastq/${sranum}" ]; then
  mkdir /scratch/ejy4bu/compBio/fastq/Old_Algae_fastq/${sranum}
fi

if [ ! -d "/scratch/ejy4bu/compBio/sra" ]; then
  mkdir /scratch/ejy4bu/compBio/sra
fi



if ls /scratch/ejy4bu/compBio/fastq/Old_Algae_fastq/${sranum}/*fastq.gz 1> /dev/null 2>&1; then
    echo "files do exist"
else
  echo "files do not exist"

  echo "force re-download"
  prefetch \
  -o /scratch/ejy4bu/compBio/sra/${sranum}.sra \
  -p \
  ${sranum}



  fasterq-dump \
  --split-files \
  --split-3 \
  -O /scratch/ejy4bu/compBio/fastq/Old_Algae_fastq/${sranum}/${sranum} \
  -e 10 \
  -p \
  --temp /scratch/ejy4bu/tmp \
  /scratch/ejy4bu/compBio/sra/${sranum}.sra

  ls -lh /scratch/ejy4bu/compBio/fastq/Old_Algae_fastq/${sranum}

fi

if [ -f "/scratch/ejy4bu/compBio/fastq/Old_Algae_fastq/${sranum}/${sranum}_1.fastq" ]; then
  gzip /scratch/ejy4bu/compBio/fastq/Old_Algae_fastq/${sranum}/${sranum}_1.fastq
  gzip /scratch/ejy4bu/compBio/fastq/Old_Algae_fastq/${sranum}/${sranum}_2.fastq
fi

if [ -f "/scratch/ejy4bu/compBio/fastq/Old_Algae_fastq/${sranum}" ]; then
  gzip -c /scratch/ejy4bu/compBio/fastq/Old_Algae_fastq/${sranum} > /scratch/ejy4bu/compBio/fastq/Old_Algae_fastq/${sranum}.fastq.gz
  rm /scratch/ejy4bu/compBio/fastq/Old_Algae_fastq/${sranum}
fi

if [ -f "/scratch/ejy4bu/compBio/sra/${sranum}.sra" ]; then 
  rm /scratch/ejy4bu/compBio/sra/${sranum}.sra
fi
#cat /home/ejy4bu/daphnia_algae/data/runs.csv | nl | grep "SRR12463313"
