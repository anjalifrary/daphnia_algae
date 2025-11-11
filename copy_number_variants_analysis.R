#ijob -A berglandlab -c2 -p standard --mem=40G
#export R_LIBS_USER="/sfs/gpfs/tardis/home/ejy4bu/R/goolf/4.5/"
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.5.0
#module load R/4.5.0
#R

library(data.table)
library(ggplot2)
library(foreach)
library(doMC)
registerDoMC(10)

out_dir <- "/scratch/ejy4bu/compBio/bam_analysis/coverage_data"

coverage_files <- list.files(out_dir, pattern = "_5000bp.csv", recursive = TRUE, full.names = TRUE)

coverage <- rbindlist(lapply(coverage_files, fread), use.names=TRUE, fill=TRUE)

coverage[, sampleID := tstrsplit(basename(V1), "_5000bp.csv")[[1]]]

coverage[, algae_source := ifelse(grepl("Old_Algae_bams", coverage_files, ignore.case = TRUE), "UTEX",
                 ifelse(grepl("Robert_samples_bams", coverage_files, ignore.case = TRUE), "REED",
                 ifelse(grepl("Sephadex", coverage_files, ignore.case = TRUE), "REED", NA)))
]