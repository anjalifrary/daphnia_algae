### in Rstudio at ood.hpc.virginia.edu 


### Documentation
# https://bioconductor.org/packages/devel/bioc/vignettes/cn.mops/inst/doc/cn.mops.pdf

# BiocManager::install("cn.mops")

library(cn.mops)
library(Rsamtools)
library(data.table)

### Output directory
out_dir <- "/scratch/ejy4bu/compBio/cnv/cnmops_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

### get input bams
sephadex_bams <- c("/scratch/ejy4bu/compBio/cnv/megabams/REED_Sephadex.megabam.bam")
REED_bams <- c("/scratch/ejy4bu/compBio/cnv/megabams/REED_NotSephadex.megabam.bam")
UTEX_bams <- c("/scratch/ejy4bu/compBio/cnv/megabams/UTEX.megabam.bam")

bam_files <- c(REED_bams, UTEX_bams)
sample_names <- c("REED_NotSephadex", "UTEX")

message("Input bams set")

### Set chlorella scaffolds
chr <- c("SIDB01000001.1","SIDB01000002.1","SIDB01000003.1","SIDB01000004.1",
         "SIDB01000005.1","SIDB01000006.1","SIDB01000007.1","SIDB01000008.1",
         "SIDB01000009.1","SIDB01000010.1","SIDB01000011.1","SIDB01000012.1",
         "SIDB01000013.1","SIDB01000014.1")

# Get read counts
bamDataRanges <- getReadCountsFromBAM(
    BAMFiles = bam_files,
    sampleNames = sample_names,
    refSeqName = chr,
    WL = 1000
)

### Run cn.mops
message("Running cn.mops")
res <- cn.mops(bamDataRanges, minWidth=2, parallel=4)
cnv_result <- calcIntegerCopyNumbers(res)

# Check what we got
cnvrs <- cnvr(cnv_result)
message("Number of CNVRs detected: ", length(cnvrs))

# Save CNV calls
cnv_df <- as.data.frame(cnvr(cnv_result))
fwrite(cnv_df, file.path(out_dir, "cnmops_CNVs.csv"), sep=",", quote=FALSE)
message("CNV calls saved. dimensions = ", nrow(cnv_df), " rows")

# Save CNVresults object
saveRDS(cnv_result, file=file.path(out_dir, "cnmops_results.rds"))

# ### Extract copy number data for plotting
# cn_data <- as.data.frame(integerCopyNumber(cnv_result))
# message("Copy number matrix: ", nrow(cn_data), " windows")
# message("Columns: ", paste(colnames(cn_data), collapse=", "))

# # Convert CN factors to numeric (e.g., "CN2" -> 2)
# cn_data$REED_NotSephadex <- as.numeric(gsub("CN", "", cn_data$REED_NotSephadex))
# cn_data$UTEX <- as.numeric(gsub("CN", "", cn_data$UTEX))

# # Add position information
# plot_data <- data.table(
#     chrom = as.character(seqnames(bamDataRanges)),
#     pos = start(bamDataRanges),
#     cn_data
# )

# # Check the data
# message("Copy number range - REED: ", min(plot_data$REED_NotSephadex), " to ", max(plot_data$REED_NotSephadex))
# message("Copy number range - UTEX: ", min(plot_data$UTEX), " to ", max(plot_data$UTEX))

# ### Plot first chromosome as example
# chr1_data <- plot_data[chrom == "SIDB01000001.1"]

# par(mfrow=c(2,1), mar=c(4,4,3,1))

# plot(chr1_data$pos, chr1_data$REED_NotSephadex, 
#      type="l", col="cyan3", lwd=2,
#      main="REED_NotSephadex - SIDB01000001.1",
#      xlab="Position (bp)", ylab="Copy Number")
# abline(h=2, col="darkgreen", lty=2, lwd=2)

# plot(chr1_data$pos, chr1_data$UTEX, 
#      type="l", col="dodgerblue3", lwd=2,
#      main="UTEX - SIDB01000001.1",
#      xlab="Position (bp)", ylab="Copy Number")
# abline(h=2, col="darkgreen", lty=2, lwd=2)

# # Summary statistics
# message("\nCNV Summary by sample:")
# for(samp in c("REED_NotSephadex", "UTEX")) {
#     gains <- sum(plot_data[[samp]] > 2)
#     losses <- sum(plot_data[[samp]] < 2)
#     normal <- sum(plot_data[[samp]] == 2)
#     message(samp, ": Gains=", gains, " Losses=", losses, " Normal=", normal)
# }

# message("done")

# # # Plot normalized data
# # plot(res, which=1)

# # # Plot segmentation
# # plot(cnv_result, which="segments")

# # # Plot individual CNVRs
# # cnvrs <- cnvr(cnv_result)
# # message("Number of CNVRs detected: ", length(cnvrs))

# # if(length(cnvrs) > 0) {
# #     # Plot first few CNVRs
# #     for(i in 1:min(5, length(cnvrs))) {
# #         plot(cnv_result, which=i)
# #         Sys.sleep(1)  # Pause between plots
# #     }
# # }

# # # # Plot 1: Raw normalized data (before CN calling)
# # # message("Creating normalized data plot...")
# # # pdf(file.path(out_dir, "cnmops_normalized.pdf"), width=14, height=8)
# # # plot(res, which=1)
# # # dev.off()

# # # # Plot 2: Segmentation results
# # # message("Creating segmentation plot...")
# # # pdf(file.path(out_dir, "cnmops_segments.pdf"), width=14, height=8)
# # # plot(cnv_result, which="segments")
# # # dev.off()

# # # # Plot 3: Individual CNVRs (if any detected)
# # # message("Plotting CNVRs...")
# # # cnvrs <- cnvr(cnv_result)
# # # if(length(cnvrs) > 0) {
# # #     # Plot first 5 CNVRs
# # #     n_plot <- min(5, length(cnvrs))
# # #     for(i in 1:n_plot) {
# # #         pdf(file.path(out_dir, paste0("cnmops_CNVR_", i, ".pdf")), width=14, height=8)
# # #         plot(cnv_result, which=i)
# # #         dev.off()
# # #         message("Plotted CNVR ", i)
# # #     }
# # # } else {
# # #     message("No CNVRs detected")
# # # }

# # # # Plot 4: Try plotting the cnvr object directly
# # # message("Plotting CNVR overview...")
# # # pdf(file.path(out_dir, "cnmops_cnvr_overview.pdf"), width=14, height=8)
# # # if(length(cnvrs) > 0) {
# # #     plot(cnvrs)
# # # }
# # # dev.off()

# # message("end of analysis")