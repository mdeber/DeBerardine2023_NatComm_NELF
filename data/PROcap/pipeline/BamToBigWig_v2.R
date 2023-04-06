#!/usr/bin/env Rscript
######################################## ########################################
## Usage: Rscript [options] <bamfile directory>
## About: This script will quality filter and modify bam files for final export
##        as bigWig files
## Dependencies include a number of Bioconductor R packages; this is  satisfied
## by installing the Bioconductor package rtracklayer
######################################## ########################################

t.start <- Sys.time() # time run-time

suppressPackageStartupMessages(
  suppressWarnings({
    library(optparse)
    library(parallel)
  })
)

msgout = function(...) {
  write(paste(...), stdout())
}

usage = "Usage: Rscript [options] -i <bamfile directory>"

option_list = list(
  make_option(c("-i", "--input-dir"),
              type = "character",
              default = getwd(),
              help = "Directory containing bam files"),
  make_option(c("-o", "--output-dir"), # output directory
              type = "character",
              default = getwd(),
              help = "Directory to output bigWig files"),
  make_option(c("-g", "--genome"),
              type = "character",
              default = "experimental",
              help = "Name to use for readcount metadata"),
  make_option(c("--metadata-output"),
              type = "character",
              default = "readcounts_BamToBigWig.txt",
              help = "File to write readcount data table"),
  make_option(c("-r", "--revcomp"), # rev comp
              type = "logical", 
              default = F, 
              action = "store_true",
              help = "Set to TRUE to take the reverse complement reads"),
  make_option(c("-e", "--read-end"), # read end
              type = "character",
              default = "neither", 
              help = "Set to '5p' or '3p' to take only that end of the read. (Applied after --revcomp)"),
  make_option(c("-s", "--stranded"), # stranded
              type = "logical",
              default = F, 
              action = "store_true",
              help = "Set to FALSE to ignore strand information"),
  make_option(c("--offset"), # offset
              type = "integer",
              default = 0, 
              help = "Offset to apply to each read. (Applied after --end and --revcomp"),
  make_option(c("--spike-genome"), # spike-in genome
              type = "character",
              default = NA,
              help = "Name of the spike-in genome to grep to find chromosome names"),
  make_option(c("--spike-out"), # spike-in output directory
              type = "character",
              default = "/dev/null/",
              help = "To output spike-in bigWig files, give an output directory"),
  make_option(c("-q", "--quality"), # quality cutoff
              type = "integer",
              default = 3, 
              help = "Minimum mapq score (the -log10 probability that the read originated from another site)"),
  make_option(c("--threads"), # threads
              type = "integer",
              default = detectCores(logical = F),
              help = "Number of threads to use for parallel processing. The default is the number of physical cores")
)

opt = parse_args(OptionParser(option_list = option_list,
                              usage = usage),
                 convert_hyphens_to_underscores = TRUE)

######################################## ########################################

### Input files ###

sample_names <- dir(opt$input_dir)
sample_names <- sample_names[grep(".bam", sample_names)]

### Ensure directories end in slashes ###

check.dirname <- function(dirname) {
  n <- nchar(dirname)
  if (substr(dirname, n, n) != "/") {
    dirname <- paste0(dirname, "/")
  }
  return(dirname)
}

opt$input_dir <- check.dirname(opt$input_dir)
opt$output_dir <- check.dirname(opt$output_dir)
opt$spike_out <- check.dirname(opt$spike_out)




msgout(" ")
msgout("########################################")
msgout("############ Bam to bigWig #############")
msgout("########################################")
msgout(" ")
msgout(" --- Input Files ---  ")
msgout(" ")
msgout(paste0("   ", sample_names), sep = "\n")
msgout(" ")
msgout(" ----- Options -----  ")
msgout(" ")
msgout(paste(paste0("   ", "--input-dir=", opt$input_dir),
             paste0("   ", "--output-dir=", opt$output_dir),
             paste0("   ", "--genome=", opt$genome),
             paste0("   ", "--metadata-output=", opt$metadata_output),
             paste0("   ", "--revcomp=", opt$revcomp),
             paste0("   ", "--read-end=", opt$read_end),
             paste0("   ", "--stranded=", opt$stranded),
             paste0("   ", "--offset=", opt$offset),
             paste0("   ", "--spike-genome=", opt$spike_genome),
             paste0("   ", "--spike-out=", opt$spike_out),
             paste0("   ", "--quality=", opt$quality),
             paste0("   ", "--threads=", opt$threads),
             sep = "\n"))
msgout(" ")
msgout(" -- Files Written --  ")
msgout(" ")

suppressPackageStartupMessages({
  suppressWarnings({
    library(Rsamtools)
    library(GenomicFiles)
    library(GenomicAlignments)
    library(rtracklayer)
    library(BiocParallel)
  })
})

######################################## ########################################

btbw <- function(bam_file, opt) {
  # Load Bam File
  infile <- BamFile(paste0(opt$input_dir, bam_file))
  param <- ScanBamParam(mapqFilter = opt$quality)
  alignment = readGAlignments(infile, 
                              use.names = F, 
                              param = param)
  
  # Initialize metadata dataframe
  readcounts.data <- data.frame(sample_name = sub(".bam", "", bam_file),
                                total_reads = NA,
                                pass_mapq = NA)
  
  # Write to metadata (in function, a 1 row dataframe)
  count.infile <- countBam(infile)
  readcounts.data$total_reads[1] <- count.infile$records
  readcounts.data$pass_mapq[1] <- length(alignment)
  
  # Make GRanges
  reads <- GRanges(alignment) 
  
  # Apply Options
  if(opt$revcomp) {
    strand(reads) = ifelse(strand(reads) == '+', '-', '+')
  }
  if(!opt$stranded ) {
    strand(reads) = "*"
  }
  if(opt$read_end == "5p" ) {
    reads = resize(reads, width = 1, fix = "start")
  } else if (opt$read_end == "3p") {
    reads = GenomicRanges::resize(reads, width = 1, fix = "end")
  }
  if(opt$offset != 0) {
    reads = GenomicRanges::shift(reads, opt$offset)
  }
  
  # Collapse identical reads
  collapse_reads <- function(reads) {
    reads <- sort(reads)
    reads.unique <- unique(reads)
    score(reads.unique) <- countOverlaps(reads.unique, reads, type = "equal")
    return(reads.unique)
  }
  reads <- collapse_reads(reads)
  
  # Base filename for bigWigs
  outfname = sub(".bam", "", bam_file, fixed = T)
  
  # If not separating spike-in reads
  if (is.na(opt$spike_genome)) {
    if(opt$stranded) {
      reads.p <- reads[which(strand(reads)=="+")]
      reads.m <- reads[which(strand(reads)=="-")]
      export.bw(reads.p, paste0(opt$output_dir, outfname, "_plus.bw"))
      export.bw(reads.m, paste0(opt$output_dir, outfname, "_minus.bw"))
      msgout("  ", paste0(opt$output_dir, outfname, "_plus.bw"))
      msgout("  ", paste0(opt$output_dir, outfname, "_minus.bw"))
    } else {
      export.bw(reads, paste0(opt$output_dir, outfname, ".bw"))
      msgout("  ", paste0(opt$output_dir, outfname, ".bw"))
    }
  } else { # If separating spike-in reads
    
    # Find spike-in chromosomes
    reads.spike <- reads
    chrom.spike <- grep(opt$spike_genome, seqlevels(reads.spike)) # get indices
    chrom.spike <- seqlevels(reads.spike)[chrom.spike] # subset
    reads.spike <- keepSeqlevels(reads.spike,
                                 chrom.spike,
                                 pruning = "tidy")
    
    # Find other chromosomes (experimental)
    reads.exp <- reads
    chrom.exp <- which( !(seqlevels(reads.exp) %in% chrom.spike) ) # get indices
    chrom.exp <- seqlevels(reads.exp)[chrom.exp]
    reads.exp <- keepSeqlevels(reads.exp,
                               chrom.exp,
                               pruning = "tidy")
    
    
    # Write to metadata: experimental & spike-in readcounts
    readcounts.data <- cbind(readcounts.data, NA, NA)
    names(readcounts.data) <- c("sample_name",
                                "total_reads",
                                "pass_mapq",
                                opt$genome,
                                opt$spike_genome)
    readcounts.data[1, 4] <- sum(score(reads.exp))
    readcounts.data[1, 5] <- sum(score(reads.spike))
    
    if(opt$stranded) {
      # Export experimental
      reads.exp.p <- reads.exp[which(strand(reads.exp) == "+")]
      reads.exp.m <- reads.exp[which(strand(reads.exp) == "-")]
      export.bw(reads.exp.p,
                paste0(opt$output_dir, outfname, "_plus.bw"))
      export.bw(reads.exp.m,
                paste0(opt$output_dir, outfname, "_minus.bw"))
      
      msgout("  ", paste0(opt$output_dir, outfname, "_plus.bw"))
      msgout("  ", paste0(opt$output_dir, outfname, "_minus.bw"))
      
      # If export spike-in (default is to /dev/null)
      if (opt$spike_out != "/dev/null/") {
        reads.spike.p <- reads.spike[which(strand(reads.spike) == "+")]
        reads.spike.m <- reads.spike[which(strand(reads.spike) == "-")]
        export.bw(reads.spike.p, 
                  paste0(opt$spike_out, outfname, "_plus.bw"))
        export.bw(reads.spike.m, 
                  paste0(opt$spike_out, outfname, "_minus.bw"))
        
        msgout("  ", paste0(opt$spike_out, outfname, "_plus.bw"))
        msgout("  ", paste0(opt$spike_out, outfname, "_minus.bw"))
      }
    } 
    
    if(!opt$stranded) { 
      export.bw(reads.exp, paste0(opt$output_dir, outfname, ".bw"))
      msgout("  ", paste0(opt$output_dir, outfname, ".bw"))
      if (opt$spike_out != "/dev/null/") {
        export.bw(reads.spike, 
                  paste0(opt$output_dir, outfname, ".bw"))
        msgout("  ", paste0(opt$spike_out, outfname, ".bw"))
      }
    }
  }
  return(readcounts.data)
}

######################################## ########################################

param.bp <- MulticoreParam(workers = opt$threads)
readcounts <- bplapply(X = sample_names, 
                       FUN = btbw, 
                       opt = opt,
                       BPPARAM = param.bp)

# combine dataframes for each sample
readcounts.data <- do.call(rbind, readcounts)

######################################## ########################################

msgout(" ")
msgout(" -- Metadata File --  ")
msgout(" ")

write.table(readcounts.data,
            opt$metadata_output,
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            append = FALSE)

msgout("  ", opt$metadata_output)
msgout(" ")
t.diff <- Sys.time() - t.start
msgout(" BamToBigWig finished in", round( as.numeric(t.diff), 1 ), units(t.diff))
msgout(" ")





