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
library(parallel)

# ncores <- 10
# cl <- makeCluster(ncores)

### Output directory
out_dir <- "/scratch/ejy4bu/compBio/cnv/cnmops_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


### get input bams
# !! first make sure bams are indexed (.bai file exists) !!
sephadex_bams <- c("/scratch/ejy4bu/compBio/cnv/megabams/REED_Sephadex.megabam.bam")
REED_bams <- c("/scratch/ejy4bu/compBio/cnv/megabams/REED_NotSephadex.megabam.bam")
UTEX_bams <- c("/scratch/ejy4bu/compBio/cnv/megabams/UTEX.megabam.bam")

bam_files <- c(sephadex_bams, REED_bams, UTEX_bams)
sample_names <- c("REED_Sephadex", "REED_NotSephadex", "UTEX")


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
#res <- cn.mops(bamDataRanges, minWidth=2)
res <- pcmops(bamDataRanges, minWidth = 2, cluster = cl)
# stopCluster(cl)


pdf(file.path(out_dir, "cnmops_rawOutput.pdf"), width=14,height=8)
plot(res)
dev.off()
message("Raw depth-read values after normalization plot saved to: ", paste0(out_dir, "cnmops_rawOutput.pdf"))

### Calculate CN values
cnv_result <- calcIntegerCopyNumbers(res)

cnv_df <- as.data.frame(cnvToGRanges(cnv_result))
cnv_df <- cnv_df[, c("sampleName", setdiff(colnames(cnv_df), "sampleName"))]
fwrite(cnv_df, file.path(out_dir, "cnmops_CNVs.csv"), sep=",", quote=FALSE)
message("CNV calls saved to ", paste0(out_dir, "cnmops_CNVs.csv"))

### Plot all samples
pdf(file.path(out_dir, "cnmops_REED_UTEX_sephadex.pdf"), width=14,height=8)
plot(cnv_result)
dev.off()

### Plot REED vs UTEX
pdf(file.path(out_dir, "cnmops_REED_UTEX.pdf"), width = 14, height = 8)
plot(cnv_result,sample=c("REED_NotSephadex", "UTEX"))
dev.off()

### Plot each group individually
for (i in seq_along(sample_names)) {
    pdf(file.path(out_dir, paste0("cnmops_", sample_names[i], ".pdf")), width=14,height=8)
    plot(cnv_result,sample=i)
    dev.off()
    message("Plot saved for ", paste0(sample_names[i]))
}
message("All plots saved")
