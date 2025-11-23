



# make conda environment with repeatmodeler in it

#conda create -n repeatmodeler_env -c bioconda -c conda-forge repeatmodeler perl=5.22.0
conda activate repeatmodeler_new

#update
wd="/scratch/rjp5nc/removedups/us_dpulex/"
cd ${wd}


#ref="/scratch/rjp5nc/Reference_genomes/orig_ref/us_pulex_ref_GCF_021134715.1_ASM2113471v1_genomic.fna"
#cleaned="/scratch/rjp5nc/removedups/us_dpulex/cleaned_US_pulex.fasta"

#seqtk seq $ref > $cleaned

#make blast db from cleaned reference genome (type nucleic acid)
#makeblastdb -in $cleaned -out my_db -dbtype nucl -title my_db -parse_seqids

#BuildDatabase -name my_db $cleaned

RepeatModeler -database my_db -pa 10