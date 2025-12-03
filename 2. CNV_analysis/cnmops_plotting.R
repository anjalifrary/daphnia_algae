library(cn.mops)
library(data.table)
library(ggplot2)
library(reshape2)

# load saved files
out_dir <- "/scratch/ejy4bu/compBio/cnv/cnmops_output"
cnv_result <- readRDS(file.path(out_dir, "cnmops_results.rds"))

# load bam data ranges
bamDataRanges <- readRDS("/scratch/ejy4bu/compBio/cnv/bamDataRanges.rds")

# get integer copy numbers
cn_data <- as.data.frame(integerCopyNumber(cnv_result))
cn_data$chrom <- as.character(seqnames(bamDataRanges))
cn_data$pos <- start(bamDataRanges)

coverage_plot <- file.path(out_dir, "cnmops_output.pdf")
  pdf(coverage_plot, width=15, height=10)
  print(
  ggplot() +
  geom_line(data=cn_data, aes(x=pos, y=REED_NotSephadex), color="cyan3") +
  geom_line(data=cn_data, aes(x=pos, y=UTEX), color="dodgerblue3") +
  geom_hline(yintercept=2, linetype="dashed", color="darkgreen") +
  facet_wrap(~chrom, scales="free_x", ncol=1) +
  labs(title="cn.mops integer copy number", x="Position (bp)", y="Copy Number") +
  theme_bw()
  )
  dev.off()
message("cnmops data saved to: ", coverage_plot)