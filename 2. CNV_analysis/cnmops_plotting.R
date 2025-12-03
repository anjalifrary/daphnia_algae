library(cn.mops)
library(data.table)
library(ggplot2)

# load saved files
out_dir <- "/scratch/ejy4bu/compBio/cnv/cnmops_output"
cnv_result <- readRDS(file.path(out_dir, "cnmops_results.rds"))

# load bam data ranges
bamDataRanges <- readRDS("/scratch/ejy4bu/compBio/cnv/cnmops_output/bamDataRanges.rds")

# get integer copy numbers
cn_data <- as.data.frame(integerCopyNumber(cnv_result))
cn_data$chrom <- as.character(seqnames(bamDataRanges))
cn_data$pos <- start(bamDataRanges)

cn_data$REED_NotSephadex <- as.numeric(gsub("CN", "", cn_data$REED_NotSephadex))
cn_data$UTEX <- as.numeric(gsub("CN", "", cn_data$UTEX))
colnames(cn_data)[colnames(cn_data) == "REED_NotSephadex"] <- "REED"


for(sc in unique(cn_data$chrom)) {
  sc_data <- cn_data[cn_data$chrom == sc, ]
  
  # File path for this scaffold
  pdf_file <- file.path(out_dir, paste0("cnmops_", sc, ".pdf"))
  pdf(pdf_file, width=20, height=10)
  
  # Plot UTEX and REED cn 
  # green line is reference (normalized at 1 for haploid chlorella)
  # each scaffold is a gray bar 

  cn_plot <- ggplot() +
    geom_line(data=sc_data, aes(x=pos, y=REED), color="cyan3", linewidth=1.2) +
    geom_line(data=sc_data, aes(x=pos, y=UTEX), color="dodgerblue3", linewidth=1.2) +
    geom_hline(yintercept=1, linetype="dashed", color="darkgreen", linewidth=1) +
    scale_y_continuous(breaks=scales::pretty_breaks(n=5)) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle=45, hjust=1, size=14),
      axis.text.y = element_text(size=14),
      axis.title.x = element_text(size=16),
      axis.title.y = element_text(size=16),
      strip.background = element_rect(fill="lightgray", color="black"),
      strip.text = element_text(face="bold", size=14)
    ) +
    labs(
      x = "Position (bp)",
      y = "Copy Number",
      title = paste("cn.mops Integer Copy Number -", sc)
    )

  print(cn_plot)
  dev.off()
  message(paste0("Saved PDF for scaffold ", sc, ": ", pdf_file))
}
