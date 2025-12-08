# set directories
in_dir <- "/scratch/ejy4bu/compBio/cnv/cnmops_output/dataFiles/haploid_5kb"
out_dir <- "/scratch/ejy4bu/compBio/cnv/cnmops_output/plots/haploid_5kb"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

# load saved cnv results
cnv_result <- readRDS(file.path(in_dir, "cnmops_results.rds"))
cnv_result <- cnvs(cnv_result)

# load bam data ranges
# bamDataRanges <- readRDS(file.path(in_dir, "bamDataRanges.rds"))

df <- as.data.frame(cnv_result)
bed <- df[,c("seqnames","start","end","CN")]

colnames(bed) <- c("chrom","start","end","name")

bed_file <- file.path(out_dir, "cnmops_calls.bed")
write.table(bed,
            file=bed_file,
            sep="\t",
            quote=FALSE,
            row.names=FALSE,
            col.names=FALSE)


message("BED file saved at:", bed_file, "\n")
