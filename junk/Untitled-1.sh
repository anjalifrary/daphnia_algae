
# Map unassembled reads
bwa mem -t 10 -K 100000000 -Y ${ref_path} \
${outfq}/${samp}/${samp}.unassembled.forward.fastq  \
${outfq}/${samp}/${samp}.unassembled.reverse.fastq | \
samtools view -Suh -q 20 -F 0x100 | \
samtools sort --threads 10 -o ${outfq}/${samp}.filt.unassembled.sort.bam
samtools index ${outfq}/${samp}.filt.unassembled.sort.bam

# Merge assembled and unassembled BAM files
samtools merge ${outfq}/${samp}.filt.merged.bam \
    ${outfq}/${samp}.sort.bam \
    ${outfq}/${samp}.filt.unassembled.sort.bam

# Index merged BAM
samtools index ${outfq}/${samp}.filt.merged.bam

# Mark duplicates
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    REMOVE_DUPLICATES=true \
    INPUT=${outfq}/chlorella_Reed.sort.bam \
    OUTPUT=${outfq}/chlorella_Reed_finalmap.bam \
    METRICS_FILE=${outfq}/${samp}_finalmap_mdups.metrics \
    CREATE_INDEX=true

# Move final BAM to output directory
mv ${outfq}/${samp}_finalmap* ${outbam}/

# Remove intermediate files
rm -f ${outfq}/${samp}.*

echo "Finished processing ${samp}"