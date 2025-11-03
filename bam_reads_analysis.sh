#!/bin/bash

BAM_DIR="/scratch/ejy4bu/compBio/bams"

CHLORELLA_LIST="/scratch/ejy4bu/compBio/genomefiles/ChrScaffoldList"

OUTFILE="/scratch/ejy4bu/compBio/bam_analysis/bam_mapping_summary.txt"
mkdir -p $OUTFILE


# Write header to output
echo -e "Sample\tTotalReads\tMappedReads\tChlorellaReads\tPctChlorella" > "$OUTPUT"

# Loop over all deduplicated BAMs
for bam in "$BAM_DIR"/*/*/*.dedup.bam; do
    sample=$(basename "$bam" .dedup.bam)
    
    # Total reads (mapped + unmapped)
    total=$(samtools idxstats "$bam" | awk '{sum+=$3+$4} END {print sum}')
    
    # Total mapped reads
    mapped=$(samtools idxstats "$bam" | awk '{sum+=$3} END {print sum}')
    
    # Reads mapped to Chlorella scaffolds
    chlorella=$(samtools idxstats "$bam" | grep -F -f "$CHLORELLA_LIST" | awk '{sum+=$3} END {print sum}')
    
    # Percentage of mapped reads assigned to Chlorella
    if [ "$mapped" -gt 0 ]; then
        pct=$(awk -v c="$chlorella" -v m="$mapped" 'BEGIN {printf "%.2f", (c/m)*100}')
    else
        pct=0
    fi
    
    # Append results to output file
    echo -e "${sample}\t${total}\t${mapped}\t${chlorella}\t${pct}" >> "$OUTFILE"
done

echo "Summary written to $OUTFILE"