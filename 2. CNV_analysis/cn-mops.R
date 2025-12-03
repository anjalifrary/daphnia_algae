# run using cn-mops.sh

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

sample_col <- intersect(colnames(cnv_df), c("sample", "sampleName"))
if(length(sample_col) > 0) {
    cnv_df <- cnv_df[, c(sample_col, setdiff(colnames(cnv_df), sample_col))]
}


## extract normalized counts
# normalized_counts <- as.data.frame(integerCopyNumber(cnv_result))
# message("Normalized count matrix dimensions: ", paste(dim(normalized_counts), collapse=" x "))



# ### Plot all samples
message("Creating plot for all samples...")
pdf(file.path(out_dir, "cnmops_REED_UTEX.pdf"), width=14,height=8)
tryCatch({
    plot(res, which=1)
    message("All samples plotted")
}, error=function(e) {
    message("Failed to plot all samples ", e$message)
})
dev.off()

### Plot each group individually
# for (s in sample_names) {
#     message("Creating plot for ", s , "...")
#     pdf(file.path(out_dir, paste0("cnmops_", s, ".pdf")), width=14, height=8)
#     tryCatch({
#         plot(cnv_result, sample=s)
#         message("Plot saved for ", s)
#     }, error = function(e){
#         message("Failed to plot for ", s, ": ", e$message)
#     })
#     dev.off()
# }
# message("All plots saved")

