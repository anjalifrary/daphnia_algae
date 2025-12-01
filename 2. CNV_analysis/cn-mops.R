#ijob -A berglandlab -c2 -p standard --mem=40G
#export R_LIBS_USER="/sfs/gpfs/tardis/home/ejy4bu/R/goolf/4.5/"
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.5.0
#module load R/4.5.0
#R

### Documentation
# https://bioconductor.org/packages/devel/bioc/vignettes/cn.mops/inst/doc/cn.mops.pdf

# BiocManager::install("cn.mops")

library(cn.mops)
library(Rsamtools)
library(data.table)

# plotting for headless environment
options(bitmapType='cairo')  # Use cairo for bitmap graphics

### Output directory
out_dir <- "/scratch/ejy4bu/compBio/cnv/cnmops_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

### get input bams
# !! first make sure bams are indexed (.bai file exists) !!
sephadex_bams <- c("/scratch/ejy4bu/compBio/cnv/megabams/REED_Sephadex.megabam.bam")
REED_bams <- c("/scratch/ejy4bu/compBio/cnv/megabams/REED_NotSephadex.megabam.bam")
UTEX_bams <- c("/scratch/ejy4bu/compBio/cnv/megabams/UTEX.megabam.bam")

bam_files <- c(REED_bams, UTEX_bams)
sample_names <- c("REED_NotSephadex", "UTEX")
#bam_files <- c(sephadex_bams, REED_bams, UTEX_bams)
#sample_names <- c("REED_Sephadex", "REED_NotSephadex", "UTEX")

message("Input bams set")

### get input bam lengths
# Set chlorella scaffolds
    chr <- c("SIDB01000001.1","SIDB01000002.1","SIDB01000003.1","SIDB01000004.1",
    "SIDB01000005.1","SIDB01000006.1","SIDB01000007.1","SIDB01000008.1","SIDB01000009.1",
    "SIDB01000010.1","SIDB01000011.1","SIDB01000012.1","SIDB01000013.1","SIDB01000014.1")

# should i specify window length? WL = ?? 1 kb windows?
bamDataRanges <- getReadCountsFromBAM(
    BAMFiles = bam_files,
    sampleNames = sample_names,
    refSeqName = chr,
    WL = 1000
)

### Run cn.mops
message("Running cn.mops")
# minWidth = minimal number of consecutive windows required for a CNV
res <- cn.mops(bamDataRanges, minWidth=2, parallel=4)
# Calculate CN values
cnv_result <- calcIntegerCopyNumbers(res)


cnv_df <- as.data.frame(cnvr(cnv_result))
fwrite(cnv_df, file.path(out_dir, "cnmops_CNVs.csv"), sep=",", quote=FALSE)
message("CNV calls saved to ", paste0(out_dir, "cnmops_CNVs.csv"))

### data for plotting
cn_matrix <- integerCopyNumber(cnv_result)
plot_data <- data.table(
    chrom = as.character(seqnames(bamDataRanges)),
    pos = start(bamDataRanges),
    REED_NotSephadex = cn_matrix[, "REED_NotSephadex"],
    UTEX = cn_matrix[, "UTEX"]
)

### plotting each chromosome
chr_list <- unique(plot_data$chrom)

for(chr_name in chr_list){
    chr_data <- plot_data[chrom == chr_name]
    
    chr_plot <- file.path(out_dir, paste0("CNV_", chr_name, ".pdf"))
    pdf(chr_plot, width=14, height=8)
    
    par(mfrow=c(2,1), mar=c(4,4,3,1))
    
    # REED plot
    plot(chr_data$pos, chr_data$REED_NotSephadex, 
         type="l", col="cyan3", lwd=2,
         main=paste("REED_NotSephadex -", chr_name),
         xlab="Position (bp)", ylab="Copy Number",
         ylim=c(0, max(chr_data$REED_NotSephadex, 4)))
    abline(h=2, col="darkgreen", lty=2, lwd=2)
    
    # UTEX plot
    plot(chr_data$pos, chr_data$UTEX, 
         type="l", col="dodgerblue3", lwd=2,
         main=paste("UTEX -", chr_name),
         xlab="Position (bp)", ylab="Copy Number",
         ylim=c(0, max(chr_data$UTEX, 4)))
    abline(h=2, col="darkgreen", lty=2, lwd=2)
    
    dev.off()
    message("Plot saved: ", chr_plot)
}

### Comparison plot
comp_plot <- file.path(out_dir, "CNV_comparison.pdf")
pdf(comp_plot, width=14, height=10)

for(chr_name in chr_list){
    chr_data <- plot_data[chrom == chr_name]
    
    plot(chr_data$pos, chr_data$REED_NotSephadex, 
         type="l", col="cyan3", lwd=2,
         main=chr_name,
         xlab="Position (bp)", ylab="Copy Number",
         ylim=c(0, max(c(chr_data$REED_NotSephadex, chr_data$UTEX), 4)))
    lines(chr_data$pos, chr_data$UTEX, col="dodgerblue3", lwd=2)
    abline(h=2, col="darkgreen", lty=2, lwd=2)
    legend("topright", legend=c("REED_NotSephadex", "UTEX"), 
           col=c("cyan3", "dodgerblue3"), lwd=2, bty="n")
}

dev.off()
message("Comparison plot saved: ", comp_plot)

message("cnmops complete")

# sample_col <- intersect(colnames(cnv_df), c("sample", "sampleName"))
# if(length(sample_col) > 0) {
#     cnv_df <- cnv_df[, c(sample_col, setdiff(colnames(cnv_df), sample_col))]
# }


## extract normalized counts
# normalized_counts <- as.data.frame(integerCopyNumber(cnv_result))
# message("Normalized count matrix dimensions: ", paste(dim(normalized_counts), collapse=" x "))



# ### Plot all samples
# message("Creating plot for all samples...")
# pdf(file.path(out_dir, "cnmops_REED_UTEX.pdf"), width=14,height=8)
# tryCatch({
#     plot(cnvr(cnv_result))
#     message("All samples plotted")
# }, error=function(e) {
#     message("Failed to plot all samples ", e$message)
# })
# dev.off()

# # ### Plot REED vs UTEX
# # pdf(file.path(out_dir, "cnmops_REED_UTEX.pdf"), width = 14, height = 8)
# # #class(res)
# # plot(cnvr(cnv_result))
# # sample=c("REED_NotSephadex", "UTEX"), useDevice=FALSE)

# #plot(res, which=1, sample = c("REED_NotSephadex", "UTEX"))
# #dev.off()

# ### Plot each group individually
# for (s in sample_names) {
#     message("Creating plot for ", s , "...")
#     pdf(file.path(out_dir, paste0("cnmops_", s, ".pdf")), width=14, height=8)
#     tryCatch({
#         plot(cnvr(cnv_result)[,s, drop=FALSE])
#         message("Plot saved for ", s)
#     }, error = function(e){
#         message("Failed to plot for ", s, ": ", e$message)
#     })
#     dev.off()
# }
# message("All plots saved")
