#ijob -A berglandlab -c2 -p standard --mem=40G
#export R_LIBS_USER="/sfs/gpfs/tardis/home/ejy4bu/R/goolf/4.5/"
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.5.0
#module load R/4.5.0
#R

library(data.table)
library(ggplot2)
library(foreach)
library(doMC)
registerDoMC(10)

out_dir <- "/scratch/ejy4bu/compBio/bam_analysis/coverage_plots"
metadata_file <- "/scratch/ejy4bu/compBio/bam_analysis/metadata.csv"

meta <- fread(metadata_file)

coverage_files <- list.files(out_dir, pattern = "_5000bp.csv", recursive = TRUE, full.names = TRUE)

coverage <- rbindlist(lapply(coverage_files, function(f) {
  dt <- fread(f)
  dt[, sampleID := sub("_5000bp\\.csv$", "", basename(f))]
  return(dt)
}), use.names = TRUE, fill = TRUE)

coverage <- merge(coverage, meta[, .(sampleID, algae_group, algae_source, propPulex)], by = "sampleID", all.x = TRUE)

setorder(coverage, algae_group, sampleID)
coverage[, sampleID := factor(sampleID, levels = unique(sampleID))]
coverage[, norm_depth := mean_depth / mean(mean_depth), by = sampleID]

### plot genome wide

genome_avg <- coverage[, .(avg_depth = mean(norm_depth)), by = .(sampleID, algae_group)]

# Dot plot
p1 <- ggplot(genome_avg, aes(x = sampleID, y = avg_depth, color = algae_group)) +
  geom_point(size = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6)) +
  labs(x = "Sample", y = "Average Depth", title = "Genome-wide mean depth by sample") +
  scale_color_manual(values = c("REED_Sephadex"="brown3",
                                "REED_NotSephadex"="cyan3",
                                "UTEX"="dodgerblue3"))

ggsave(file.path(out_dir, "genome_wide_coverage.pdf"), p1, width = 20, height = 10)


### per chromosome plots
chr_list <- unique(coverage$chrom)

for(chr in chr_list){
  chr_data <- coverage[chrom == chr]
  
  # Create numeric x-axis along the chromosome (window index)
  setorder(chr_data, start)
  chr_data[, window_mid := (start + end)/2]    

  chr_plot <- file.path(out_dir, paste0("Coverage_", chr, ".pdf"))
  pdf(chr_plot, width=20, height=10)
  print(ggplot(chr_data, aes(x = window_mid, y = norm_depth, color = algae_group, group = sampleID)) +
    geom_line() +
    #geom_smooth(se=FALSE, span=0.1)+
    theme_bw() +
    labs(x = "Genome position", 
            y = "Normalized Coverage",
            title = paste0("Normalized Coverage across ", chr)) +
    scale_color_manual(values = c("REED_Sephadex"="brown3",
                                  "REED_NotSephadex"="cyan3",
                                  "UTEX"="dodgerblue3"))
  )
  dev.off()
  message("Chromosomal coverage plot written to: ", chr_plot)
}
