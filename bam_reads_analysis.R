#ijob -A berglandlab -c2 -p standard --mem=40G
#export R_LIBS_USER="/sfs/gpfs/tardis/home/ejy4bu/R/goolf/4.5/"
#module load gcc/11.4.0  openmpi/4.1.4 icu R/4.5.0
#module load R/4.5.0
#R

### install a new package; you only need to do this once.
  #.libPaths("/sfs/gpfs/tardis/home/ejy4bu/R/goolf/4.5/")
  #cat("R library paths:", .libPaths(), "\n")
  
  if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

  #BiocManager::install("Rsamtools")

### libraries
  library(ggplot2)
  library(data.table)
  library(foreach)
  library(doMC)
  registerDoMC(10)
  library(Rsamtools)

### changed cores from 2 to 10

### this
### specify the bam file
  bam_files <- system("ls -d /scratch/ejy4bu/compBio/bams/*/*/*.dedup.bam", intern=T)
  # bam_files <- bam_files[1:10] #test a subset
  message("Found ", length(bam_files), " BAM files")

  # /scratch/ejy4bu/compBio/bams/Robert_samples_bams/RobertUK_G9/RobertUK_G9.dedup.bam"


### Set chlorella scaffolds
    chr <- c("SIDB01000001.1","SIDB01000002.1","SIDB01000003.1","SIDB01000004.1",
    "SIDB01000005.1","SIDB01000006.1","SIDB01000007.1","SIDB01000008.1","SIDB01000009.1",
    "SIDB01000010.1","SIDB01000011.1","SIDB01000012.1","SIDB01000013.1","SIDB01000014.1")

### How many reads map to each chromosome?
    reads <- foreach(bamFile=bam_files, .combine="rbind")%dopar%{
      message("Processing: ", bamFile)

    ### get the information about number of reads for each chromosome

    # ensure file exists
      if(!file.exists(paste(bamFile, ".bai", sep=""))) indexBam(bamFile)

      # get total read counts
      stats <- as.data.table(idxstatsBam(bamFile))
      setkey(stats, seqnames)


    ### subset to chlorella scaffolds
      stats_pulex <- stats[J(chr)]
      #stats_pulex[,species:="pulex"]

    ### get proportion of reads that mapped to chlorella
      prop_pulex <- sum(stats_pulex$mapped)/sum(stats$mapped)

    ### get sample ID
      sample_id <- sub(".dedup.bam", "", basename(bamFile))

    ### print data table
      data.table(sampleID = sample_id, propPulex=prop_pulex)
  }

  out_dir <- "/scratch/ejy4bu/compBio/bam_analysis"
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


  ### create metadata csv table
  out_metadata <- file.path(out_dir, "metadata.csv")

  meta <- data.table(
    sampleID = reads$sampleID,
    bam_path = bam_files,
    algae_source = ifelse(grepl("Old_Algae_bams", bam_files, ignore.case = TRUE), "UTEX",
                  ifelse(grepl("Robert_samples_bams", bam_files, ignore.case = TRUE), "REED",
                  ifelse(grepl("Sephadex", bam_files, ignore.case = TRUE), "REED", NA))),
    sephadex = ifelse(grepl("Old_Algae_bams", bam_files, ignore.case = TRUE), "N",
                  ifelse(grepl("Robert_samples_bams", bam_files, ignore.case = TRUE), "N",
                  ifelse(grepl("Sephadex", bam_files, ignore.case = TRUE), "Y", NA))),
    algae_group = ifelse(algae_source=="REED" & sephadex=="N", "REED_NotSephadex",
                  ifelse(algae_source=="REED" & sephadex=="Y", "REED_Sephadex",
                  ifelse(algae_source=="UTEX", "UTEX")))
    propPulex = reads$propPulex *100
  )

#convert algae_source to factor for plotting by color
  meta[, algae_source := factor(algae_source, levels = c("REED", "UTEX"))]

  fwrite(meta, out_metadata)
  message("Metadata written to: ", out_metadata)

  meta[, sampleID := factor(sampleID, levels = meta[order(-propPulex)]$sampleID)]

  ### write results to csv
  out_file <- file.path(out_dir, "bam_pulex_reads.csv")

  fwrite(reads, out_file)
  message("wrote results to: ", out_file)


  ### summary of results
  summary(reads$propPulex)

  ### plot data
  plot_faceted <- "/scratch/ejy4bu/compBio/bam_analysis/bam_pulex_faceted_plot.pdf"
  pdf(plot_faceted, width=12, height=6)
  print(
  ggplot(meta,
      aes(x=sampleID, y = propPulex, color = algae_group)) + 
      geom_point() + 
      facet_wrap(~algae_source) +
      ggtitle("Chlorella proportion by algae source") +
      ylab("%Chlorella") + 
      xlab("Sample ID") +
      #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_color_manual(values = c("REED_Sephadex" = "red", "REED_NotSephadex" = "cyan", "UTEX" = "blue"))
  )
  dev.off()

  ### box plot
  plot_box <- "/scratch/ejy4bu/compBio/bam_analysis/bam_pulex_box_plot.pdf"
  pdf(plot_box, width=12, height=6)
  print(
    ggplot(meta,
    aes(x = algae_group, y = propPulex, fill = algae_group)) +
    geom_boxplot(alpha = 0.7, outlier.color = "black") +
    #geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
    facet_wrap(~ algae_source, scales = "free_x") +
    ylab("% Chlorella") +
    xlab("Algae Source") +
    ggtitle("Distribution of % Chlorella by Algae Source") +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = c("REED_Sephadex" = "red", "REED_NotSephadex" = "cyan", "UTEX" = "blue"))
  )
  dev.off()




  #rd[,realProp:=nSim/(nSim+nMel)]
  #rd
  #summary(rd$propPulex)
  #fwrite(rd, "/scratch/ejy4bu/compBio/bams/bam_reads_propPulex.csv")


### DISCARDS

### estimated vs real
  #ggplot(data=rd, aes(x=realProp, y=propSim)) + geom_point()

  
      #out[,nPulex:=as.numeric(tstrsplit(bamFile, "\\.")[[8]])]
      #out[,nMel:=as.numeric(tstrsplit(bamFile, "\\.")[[10]])]
      #out[,samp:=last(tstrsplit(bamFile, "/"))]

            #stats_small[!grepl("sim", seqnames),species:="mel"]

        ### this command aggregates the data. For each unique value of species, we calculate the total number of mapped reads.
      # stats_small.ag <- stats_small[,list(mapped=sum(mapped), seqlength=sum(seqlength)), list(species)]