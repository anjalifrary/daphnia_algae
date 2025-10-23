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
  bam_files <- bam_files[1:10] #test a subset
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

  ### write results to csv
  out_dir <- "/scratch/ejy4bu/compBio/bam_analysis"
  dir.create(out_dir)

  out_file <- file.path(out_dir, "bam_pulex_reads.csv")

  fwrite(reads, out_file)
  message("wrote results to: ", out_file)

  ### summary of results
  summary(reads$propPulex)

  ### plot data
  plot_file <- "/scratch/ejy4bu/compBio/bam_analysis/bam_pulex_plot.pdf"
  pdf(plot_file, width=12, height=6)
  print(
  ggplot(reads, aes(x=sampleID, y = propPulex * 100)) + geom_point() + ylab("%Chlorella") + xlab("Sample ID")
  )

  dev.off()

  message("plotted data at: ", "/scratch/ejy4bu/compBio/bam_analysis/bam_pulex_plot.pdf")


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