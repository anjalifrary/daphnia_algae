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
#megabams <- system("ls -d /scratch/ejy4bu/compBio/cnv/megabams/*.megabam.bam", intern=T)
megabam <- "/scratch/ejy4bu/compBio/cnv/megabams/*.megabam.bam"

### Set chlorella scaffolds
    chr <- c("SIDB01000001.1","SIDB01000002.1","SIDB01000003.1","SIDB01000004.1",
    "SIDB01000005.1","SIDB01000006.1","SIDB01000007.1","SIDB01000008.1","SIDB01000009.1",
    "SIDB01000010.1","SIDB01000011.1","SIDB01000012.1","SIDB01000013.1","SIDB01000014.1")

bamfile <- BamFile(megabam)
open(bamfile)

out_dir <- "/scratch/ejy4bu/compBio/cnv/megabamCoverageAnalysis"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)



### per chromosome plots
chr_list <- unique(coverage$chrom)

for(chr in chr_list){
  chr_data <- coverage[chrom == chr]
  
  # Create numeric x-axis along the chromosome (window index)
  setorder(chr_data, start)
  chr_data[, window_mid := (start + end)/2]  

  ### calculate group mean coverage
#   chr_group_data <- chr_data[, .(
#     mean_norm_depth = mean(norm_depth, na.rm=TRUE)
#     ), by = .(window_mid, algae_group)]

  #chr_data[norm_depth > 5, norm_depth := 5]  
  #chr_group_data[mean_norm_depth > 5, mean_norm_depth := 5]

  #dir.create(file.path(out_dir), showWarnings=FALSE,recursive=TRUE)

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
