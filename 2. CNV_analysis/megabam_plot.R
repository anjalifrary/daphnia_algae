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
megabams <- Sys.glob("/scratch/ejy4bu/compBio/cnv/megabams/*.megabam.bam")

### Create metadata for groups
# Modify this depending on your filenames
sample_table <- data.table(
  bam = megabams,
  sample = basename(megabams),
  algae_group = fifelse(grepl("REED_Sephadex", megabams), "REED_Sephadex",
                 fifelse(grepl("REED_NotSephadex", megabams), "REED_NotSephadex",
                 "UTEX"))
)

### Chlorella scaffolds
scaffolds <- c("SIDB01000001.1","SIDB01000002.1","SIDB01000003.1","SIDB01000004.1",
               "SIDB01000005.1","SIDB01000006.1","SIDB01000007.1","SIDB01000008.1",
               "SIDB01000009.1","SIDB01000010.1","SIDB01000011.1","SIDB01000012.1",
               "SIDB01000013.1","SIDB01000014.1")

### Output directory
out_dir <- "/scratch/ejy4bu/compBio/cnv/megabamCoverageAnalysis"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)



### 1. FUNCTION: compute windowed coverage per bam     ###

compute_cov <- function(bam, chr, win=5000) {

  # Find chromosome length from BAM header
  hdr <- scanBamHeader(bam)[[1]]$targets
  chr_len <- hdr[[chr]]
  if (is.null(chr_len)) {
    stop(paste("Chromosome", chr, "not found in", bam))
  }

  windows <- tileGenome(seqlengths = chr_len,
                        tilewidth = win,
                        cut.last.tile.in.chrom = TRUE)

  param <- ScanBamParam(
    which = windows,
    what = c("pos")
  )

  # Per-window coverage
  cov <- coverage(BamFile(bam), param=param)

  # Extract numeric vector coverage
  chr_cov <- as.numeric(cov[[chr]])

  # convert to data.table
  dt <- data.table(
    chrom = chr,
    start = start(windows),
    end = end(windows),
    depth = chr_cov
  )

  dt[, bam := bam]

  return(dt)
}


### 2. RUN coverage on all bams for all chromosomes   ###

message("Computing coverage...\n")

coverage_dt <- rbindlist(
  foreach(b = megabams, .combine="rbind") %dopar% {
    foreach(c = scaffolds, .combine="rbind") %dopar% {
      compute_cov(b, c, win=5000)
    }
  }
)

message("Finished coverage.\n")

### Merge sample group information
coverage_dt <- merge(coverage_dt, sample_table, by="bam")



### 3. Normalization (per-sample median normalization) ###

coverage_dt[, norm_depth := depth / median(depth, na.rm=TRUE), by = bam]


### 4. Compute per-group mean coverage for plotting    ###

group_cov <- coverage_dt[, .(
  mean_norm_depth = mean(norm_depth, na.rm=TRUE)
), by = .(chrom, start, end, algae_group)]

group_cov[, window_mid := (start + end) / 2]


### 5. Plot each chromosome                           ###

for (chr in scaffolds) {

  chr_data <- group_cov[chrom == chr]

  plot_file <- file.path(out_dir, paste0("Coverage_", chr, ".pdf"))
  pdf(plot_file, width=20, height=10)

  print(
    ggplot(chr_data,
           aes(x = window_mid,
               y = mean_norm_depth,
               color = algae_group,
               group = algae_group)) +
      geom_line() +
      theme_bw() +
      labs(
        x = "Genome position",
        y = "Normalized megabam coverage",
        title = paste0("Megabam Mean Coverage: ", chr)
      ) +
      scale_color_manual(values = c(
        "REED_Sephadex"="brown3",
        "REED_NotSephadex"="cyan3",
        "UTEX"="dodgerblue3"
      ))
  )

  dev.off()
  message("Wrote: ", plot_file)
}

message("All done!\n")
