#ijob -A berglandlab -c2 -p standard --mem=40G
#export R_LIBS_USER="/sfs/gpfs/tardis/home/ejy4bu/R/goolf/4.5/"
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.5.0
#module load R/4.5.0
#R

### libraries
library(Rsamtools)
library(ggplot2)


library(data.table)
library(foreach)
library(doMC)
registerDoMC(10)


### Define megabam paths
megabams <- system("ls -d /scratch/ejy4bu/compBio/cnv/megabams/*.megabam.bam", intern=T)
#megabam <- "/scratch/ejy4bu/compBio/cnv/megabams/*.megabam.bam"
message("Found ", length(megabams), " megabam files")

out_dir <- "/scratch/ejy4bu/compBio/cnv/megabamCoverageAnalysis"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


### Set chlorella scaffolds
    chr <- c("SIDB01000001.1","SIDB01000002.1","SIDB01000003.1","SIDB01000004.1",
    "SIDB01000005.1","SIDB01000006.1","SIDB01000007.1","SIDB01000008.1","SIDB01000009.1",
    "SIDB01000010.1","SIDB01000011.1","SIDB01000012.1","SIDB01000013.1","SIDB01000014.1")


# get chromosome lengths from scaffold file
chrom_lengths <- fread("/scratch/ejy4bu/compBio/genomefiles/scaffold_lengths.txt",
                       col.names = c("chr_names", "chr_lengths"))
setkey(chrom_lengths, chr_names)

# get read length
bam <- scanBam(megabams[1])
read_length <- mean(nchar(as.character(bam[[1]]$seq)))

# calculate coverage per megaBAM per chromosome
coverage <- foreach(bamFile = megabams, .combine="rbind") %dopar% {
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
    stats_pulex[, sampleID := sub(".megabam.bam", "", basename(bamFile))]

    stats_pulex[, .(sampleID, chr_names, mapped, chr_lengths, coverage)]
}



# data table
out_file <- file.path(out_dir, "megabam_coverage_analysis.pdf")

### per chromosome plots
chr_list <- unique(coverage$chrom)

for(chr in chr_list){
  chr_data <- coverage[chrom == chr]
  
  # Create numeric x-axis along the chromosome (window index)
  setorder(chr_data, start)
  chr_data[, window_mid := (start + end)/2]  

  ### calculate group mean coverage
  chr_group_data <- chr_data[, .(
    mean_norm_depth = mean(norm_depth, na.rm=TRUE)
    ), by = .(window_mid, algae_group)]

  #chr_data[norm_depth > 5, norm_depth := 5]  
  #chr_group_data[mean_norm_depth > 5, mean_norm_depth := 5]


  chr_plot <- file.path(out_dir, paste0("Coverage_", chr, ".pdf"))
  pdf(chr_plot, width=20, height=10)
  print(ggplot(chr_group_data, aes(x = window_mid, y = mean_norm_depth, color = algae_group, group = algae_group)) +
    geom_line() +
    #geom_smooth(se=FALSE, span=0.1)+
    theme_bw() +
    labs(x = "Genome position", 
            y = "Megabam coverage",
            title = paste0("Megabam Aggregate Coverage across ", chr)) +
    scale_color_manual(values = c("REED_Sephadex"="brown3",
                                  "REED_NotSephadex"="cyan3",
                                  "UTEX"="dodgerblue3"))
  )
  dev.off()
  message("Chromosomal megabam coverage plot written to: ", chr_plot)
}
