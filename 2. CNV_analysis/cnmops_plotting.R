library(cn.mops)
library(data.table)
library(ggplot2)
library(reshape2)

# load saved files
out_dir <- "/scratch/ejy4bu/compBio/cnv/cnmops_output"
cnv_result <- readRDS(file.path(out_dir, "cnmops_results.rds"))

# 