#ijob -A berglandlab -c2 -p standard --mem=40G
#export R_LIBS_USER="/sfs/gpfs/tardis/home/ejy4bu/R/goolf/4.5/"
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.5.0
#module load R/4.5.0
#R

library(Rsamtools)
library(data.table)
library(ggplot2)
library(foreach)
library(doMC)
registerDoMC(10)



# get bam files
bam_files <- system("ls -d /scratch/ejy4bu/compBio/bams/*/*/*.dedup.bam", intern=T)
  # bam_files <- bam_files[1:10] #test a subset
  message("Found ", length(bam_files), " BAM files")

### Set chlorella scaffolds
chr <- c("SIDB01000001.1","SIDB01000002.1","SIDB01000003.1","SIDB01000004.1",
    "SIDB01000005.1","SIDB01000006.1","SIDB01000007.1","SIDB01000008.1","SIDB01000009.1",
    "SIDB01000010.1","SIDB01000011.1","SIDB01000012.1","SIDB01000013.1","SIDB01000014.1")

# get chromosome lengths from scaffold file
chrom_lengths <- fread("/scratch/ejy4bu/compBio/genomefiles/scaffold_lengths.txt",
                       col.names = c("chr_names", "chr_lengths"))
setkey(chrom_lengths, chr_names)

# get read length
bam <- scanBam(bam_files[1])
read_length <- mean(nchar(as.character(bam[[1]]$seq)))

# calculate coverage per BAM / per chromosome
coverage <- foreach(bamFile = bam_files, .combine="rbind") %dopar% {
  message("Processing: ", bamFile)
  
  if(!file.exists(paste(bamFile, ".bai", sep=""))) indexBam(bamFile)

  # get mapped reads per chromosome 
  stats <- as.data.table(idxstatsBam(bamFile))
  setnames(stats, "seqnames", "chr_names")  # rename
  setkey(stats, chr_names)

  
  # subset to chlorella chromosomes 
  stats_pulex <- stats[J(chr)]
  
  # merge with chromosome lengths
  stats_pulex <- merge(stats_pulex, chrom_lengths, by="chr_names")

  # calculate coverage: (mapped reads * read length) / chromosome length
  stats_pulex[, coverage := (mapped * read_length) / chr_lengths]

  # add sample ID
  stats_pulex[, sampleID := sub(".dedup.bam", "", basename(bamFile))]

  stats_pulex[, .(sampleID, chr_names, mapped, chr_lengths, coverage)]
}

out_dir <- "/scratch/ejy4bu/compBio/bam_analysis"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# save data table
out_file <- file.path(out_dir, "bam_coverage_table.csv")
fwrite(coverage, out_file)
message("Coverage table written to: ", out_file)

# generate coverage plot 
# coverage[, coverage := as.numeric(coverage)]
coverage_avg <- coverage[, .(coverage=mean(coverage)), by = chr_names]
coverage_avg[, chr_names := factor(chr_names, levels=chr)] # order scaffold names

coverage_avg <- file.path(out_dir, "avg_coverage_per_chromosome.pdf")
  pdf(coverage_avg, width=40, height=10)
  print(
  ggplot(coverage_avg, aes(x=chr_names, y = coverage)) + 
      geom_bar(stat="identity") +
      ylab("Average Coverage") + 
      xlab("Scaffold") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size=12))
)
  dev.off()
message("Coverage plot written to: ", coverage_avg)
