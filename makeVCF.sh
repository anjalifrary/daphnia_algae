module load samtools varscan bcftools
# sbatch --array=1-12 ~/Genoseq_2023/run_varscan_copy.sh
# 40646429
# Get the chromosome for the current SLURM task
chr=$(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/ejy4bu/compBio/genomefiles/ChrScaffoldList)
echo $chr
# Store BAM file list in a variable to maintain order
bam_list=$(ls /scratch/rjp5nc/UK2022_2024/allshortreads/sortedbamsdedup/*.bam | sort)
# Generate VCF with VarScan
samtools mpileup \
    -r ${chr} \
    --fasta-ref /project/berglandlab/daphnia_ref/totalHiCwithallbestgapclosed.fa \
    $bam_list | \
java -jar $EBROOTVARSCAN/VarScan.v2.4.4.jar mpileup2snp \
    /dev/stdin \
    --min-coverage 4 \
    --min-var-freq 0.001 \
    --output-vcf > /scratch/rjp5nc/UK2022_2024/allshortreads/chr/2022seq.${chr}.vcf
