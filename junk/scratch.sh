infq="/scratch/ejy4bu/compBio/fastq/Old_Algae_fastq"
samples=$(ls -d ${infq}/*) #Array to folder paths


ls $infq/*/_trimmedmerged1.fq.gz

for samp in $samples; do
    if ls ${samp}/_trimmedmerged1.fq.gz 1> /dev/null 2>&1; then
        rm $samp/_trimmedmerged1.fq.gz
        echo "deleted $samp/_trimmedmerged1.fq.gz"
    fi
        if ls ${samp}/_trimmedmerged2.fq.gz 1> /dev/null 2>&1; then
        rm $samp/_trimmedmerged2.fq.gz
        echo "deleted $samp/_trimmedmerged2.fq.gz"
    fi
done