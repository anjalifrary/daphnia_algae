library(cn.mops)
library(data.table)
library(ggplot2)

# load saved files
out_dir <- "/scratch/ejy4bu/compBio/cnv/cnmops_output"
cnv_result <- readRDS(file.path(out_dir, "cnmops_results.rds"))

# load bam data ranges
bamDataRanges <- readRDS("/scratch/ejy4bu/compBio/cnv/bamDataRanges.rds")

# get integer copy numbers
cn_data <- as.data.frame(integerCopyNumber(cnv_result))
cn_data$chrom <- as.character(seqnames(bamDataRanges))
cn_data$pos <- start(bamDataRanges)

cn_data$REED_NotSephadex <- as.numeric(gsub("CN", "", cn_data$REED_NotSephadex))
cn_data$UTEX <- as.numeric(gsub("CN", "", cn_data$UTEX))

coverage_plot <- file.path(out_dir, "cnmops_genomewide_plot.pdf")
pdf(coverage_plot, width=40, height=40)

ggplot() +
  geom_line(data=cn_data, aes(x=pos, y=REED_NotSephadex, color="REED_NotSephadex"), linewidth=1.2) +
  geom_line(data=cn_data, aes(x=pos, y=UTEX, color="UTEX"), linewidth=1.2) +
  geom_hline(yintercept=2, linetype="dashed", color="darkgreen", linewidth=1) +
  facet_wrap(~chrom, scales="free_x", ncol=1) +
  scale_color_manual(values=c("REED_NotSephadex"="cyan3", "UTEX"="dodgerblue3")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=16),
    axis.title.y = element_text(size=16),
    legend.position = "top",
    strip.background = element_rect(fill="lightgray", color="black"),
    strip.text = element_text(face="bold", size=14)
  ) +
  labs(
    x = "Position (bp)",
    y = "Copy Number",
    color = "Sample",
    title = "cn.mops Integer Copy Number Across Genome"
  )

dev.off()
message("saved to ", coverage_plot)