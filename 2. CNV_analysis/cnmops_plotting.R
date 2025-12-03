# run using cnmops_plotting_runRscript.sh

library(cn.mops)
library(data.table)
library(ggplot2)

# set directories
in_dir <- "/scratch/ejy4bu/compBio/cnv/cnmops_output/dataFiles/haploid_5kb"
out_dir <- "/scratch/ejy4bu/compBio/cnv/cnmops_output/plots/haploid_5kb"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# load saved cnv results
cnv_result <- readRDS(file.path(in_dir, "cnmops_results.rds"))

# load bam data ranges
bamDataRanges <- readRDS(file.path(in_dir, "bamDataRanges.rds"))

# get integer copy numbers
cn_data <- as.data.frame(integerCopyNumber(cnv_result))
cn_data$chrom <- as.character(seqnames(bamDataRanges))
cn_data$pos <- start(bamDataRanges)

cn_data$REED_NotSephadex <- as.numeric(gsub("CN", "", cn_data$REED_NotSephadex))
cn_data$UTEX <- as.numeric(gsub("CN", "", cn_data$UTEX))
colnames(cn_data)[colnames(cn_data) == "REED_NotSephadex"] <- "REED"

# continuous cn.mops result
# this function is exerimental and only works on a reference ?? 
# cn_cont <- as.data.frame(calcFractionalCopyNumbers(cnv_result))
# cn_cont$chrom <- as.character(seqnames(bamDataRanges))
# cn_cont$pos <- start(bamDataRanges)

wrap_size <- 1000  # adjust based on scaffold length

for(sc in unique(cn_data$chrom)) {
  sc_data <- cn_data[cn_data$chrom == sc, ]
  
  # each line group  for wrapping
  sc_data$line <- ceiling(seq_len(nrow(sc_data)) / wrap_size)
  
  # File path for this scaffold
  pdf_file <- file.path(out_dir, paste0("cnmops_", sc, ".pdf"))
  pdf(pdf_file, width=20, height=10)
  
  cn_plot <- ggplot(sc_data, aes(x=pos)) +
    geom_line(aes(y=REED, color="REED"), linewidth=1.2) +
    geom_line(aes(y=UTEX, color="UTEX"), linewidth=1.2) +
    facet_wrap(~line, scales="free_x", ncol=1) +  # each line is a section of the scaffold
    scale_color_manual(values=c("REED"="cyan3", "UTEX"="dodgerblue3")) +
    scale_y_continuous(breaks=scales::pretty_breaks(n=5)) +
    theme_bw() +
    theme(
        plot.title = element_text(size=22, face="bold"),
      axis.text.x = element_text(angle=45, hjust=1, size=14),
      axis.text.y = element_text(size=14),
      axis.title.x = element_text(size=16),
      axis.title.y = element_text(size=16),
      legend.position="top",
      legend.title = element_text(size=18, face="bold"),
      legend.text = element_text(size=16),
      strip.background = element_rect(fill="lightgray", color="black"),
      strip.text = element_blank()
    ) +
    labs(
      x = "Position (bp)",
      y = "Copy Number",
      color = "Sample",
      title = paste("cn.mops Integer Copy Number -", sc)
    )
  
  print(cn_plot)
  dev.off()
  message(paste0("Saved wrapped PDF for scaffold ", sc, ": ", pdf_file))
}