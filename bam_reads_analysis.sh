OUTFILE=/scratch/ejy4bu/compBio/bam_analysis/bam_mapping_summary.txt
mkdir -p $OUTFILE

for bam in /scratch/ejy4bu/compBio/bams/*/*/*.dedup.bam; do
    sample=$(basename "$bam" .dedup.bam)
    
    # Total reads
    total=$(samtools idxstats "$bam" | awk '{sum+=$3+$4} END {print sum}')
    
    # Total mapped reads
    mapped=$(samtools idxstats "$bam" | awk '{sum+=$3} END {print sum}')
    
    # Mapped reads on Chlorella scaffolds
    chlorella=$(samtools idxstats "$bam" | grep -F -f chr_list.txt | awk '{sum+=$3} END {print sum}')
    
    # % mapped to Chlorella
    pct=$(awk -v c="$chlorella" -v m="$mapped" 'BEGIN {if(m>0) print (c/m)*100; else print 0}')
    
    # Print results
    echo -e "${sample}\t${total}\t${mapped}\t${chlorella}\t${pct}"
done > $OUTFILE
