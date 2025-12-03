#ijob -A berglandlab -c4 -p standard --mem=100G
#export R_LIBS_USER="/sfs/gpfs/tardis/home/ejy4bu/R/goolf/4.5/"
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.5.0
#module load R/4.5.0
#R

library(Rsamtools)
library(data.table)
library(ggplot2)
library(foreach)
library(doMC)
registerDoMC(2)



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

bam_chunks <- split(bam_files, ceiling(seq_along(bam_files)/50))
coverage <- rbindlist(lapply(bam_chunks, function(bams){
  foreach(bamFile = bams, .combine="rbind") %dopar% {
# coverage <- foreach(bamFile = bam_files, .combine="rbind") %dopar% {
  message("Processing: ", bamFile)
  
  if(!file.exists(paste0(bamFile, ".bai", sep=""))) indexBam(bamFile)

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
}))

meta <- data.table(
  sampleID = sub(".dedup.bam", "", basename(bam_files)),
  bam_path = bam_files,
  algae_source = ifelse(grepl("Old_Algae_bams", bam_files, ignore.case = TRUE), "UTEX",
                 ifelse(grepl("Robert_samples_bams", bam_files, ignore.case = TRUE), "REED",
                 ifelse(grepl("Sephadex", bam_files, ignore.case = TRUE), "REED", NA))),
  sephadex = ifelse(grepl("Old_Algae_bams", bam_files, ignore.case = TRUE), "N",
             ifelse(grepl("Robert_samples_bams", bam_files, ignore.case = TRUE), "N",
             ifelse(grepl("Sephadex", bam_files, ignore.case = TRUE), "Y", NA)))
)

#convert algae_source to factor for plotting by color

meta[, algae_group := ifelse(algae_source=="REED" & sephadex=="N", "REED_NotSephadex",
                      ifelse(algae_source=="REED" & sephadex=="Y", "REED_Sephadex",
                      ifelse(algae_source=="UTEX", "UTEX", NA)))]

meta[, algae_source := factor(algae_source, levels = c("REED", "UTEX"))]
meta[, algae_group := factor(algae_group, levels = c("REED_NotSephadex", "REED_Sephadex", "UTEX"))]

coverage <- merge(coverage, meta[, .(sampleID, algae_group)], by = "sampleID", all.x = TRUE)
coverage[, coverage := as.numeric(coverage)]


out_dir <- "/scratch/ejy4bu/compBio/bam_analysis/coverage_plots/mean_depth"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# save data table
# out_file <- file.path(out_dir, "bam_coverage_table.csv")
# fwrite(coverage, out_file)
# message("Coverage table written to: ", out_file)



### generate coverage plot by chromosomes for two experimental groups


coverage_avg <- coverage[, .(coverage=mean(coverage)), by = chr_names]
coverage_avg[, chr_names := factor(chr_names, levels=chr)] # order scaffold names

coverage_sub <- coverage[algae_group %in% c("REED_NotSephadex", "UTEX")]

coverage_avg_sub <- coverage_sub[, .(coverage = mean(coverage)), by = .(chr_names, algae_group)]
coverage_avg_sub[, chr_names := factor(chr_names, level=chr)]

coverage_plot <- file.path(out_dir, "avg_coverage_per_chromosome_noSeph.pdf")
  pdf(coverage_plot, width=20, height=10)
  print(
  ggplot(coverage_avg,
      aes(x=chr_names, y = coverage, fill=algae_group)) + 
      geom_bar(stat="identity", position = position_dodge(width = 0.7), width = 0.3) +
      ylab("Average Coverage") + 
      xlab("Scaffold") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size=14)) + 
      scale_fill_manual(values = c("REED_NotSephadex" = "cyan3",
                                 "UTEX" = "dodgerblue3"))
)
  dev.off()
message("Grouped chr coverage plot written to: ", coverage_plot)



# ### generate coverage plot for genome wide depth by sample


# genome_avg <- coverage[, .(avg_coverage = mean(coverage)), by = .(algae_group, sampleID)]
# setorder(genome_avg, algae_group, sampleID)

# genome_avg[, sampleID := factor(sampleID, levels = genome_avg[order(algae_group, sampleID)]$sampleID)]

# genome_plot <- file.path(out_dir, "avg_coverage_across_genome.pdf")
# pdf(genome_plot, width=20, height=10)
# print(
#   ggplot(genome_avg, aes(x=sampleID, y=avg_coverage, fill = algae_group)) + 
#   ggtitle("Mean depth across genome by sample") +
#   geom_bar(stat="identity") + 
#   ylab("Mean Coverage") + 
#   xlab("Sample Name") + 
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
#   scale_fill_manual(values = c("REED_Sephadex" = "brown3", "REED_NotSephadex" = "cyan3", "UTEX" = "dodgerblue3"))
# )
# dev.off()
# message("Genome wide coverage plot written to: ", genome_plot)



# ### generate coverage plot for each sample's mean depth

# for (chrom in chr) {
#   chr_data <- coverage[chr_names == chrom]

#   chr_avg <- chr_data[, .(avg_coverage = mean(coverage)), by = .(algae_group, sampleID)]
#   setorder(chr_avg, algae_group, sampleID)

#   chr_avg[, sampleID := factor(sampleID, levels = chr_avg[order(algae_group, sampleID)]$sampleID)]  

#   chr_plot <- file.path(out_dir, paste0("coverage_", chrom, ".pdf"))
#   pdf(chr_plot, width=20, height=10)
#   print(
#     ggplot(chr_avg, aes(x=sampleID, y=avg_coverage, fill=algae_group)) + 
#     ggtitle(paste("Mean depth across", chrom, "by sample")) +
#     geom_bar(stat="identity") + 
#     ylab("Mean Coverage") + 
#     xlab("Sample Name") + 
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
#     scale_fill_manual(values = c("REED_Sephadex" = "brown3", "REED_NotSephadex" = "cyan3", "UTEX" = "dodgerblue3"))
#   )
#   dev.off()
#   message("Chromosomal coverage plot written to: ", chr_plot)
# } 