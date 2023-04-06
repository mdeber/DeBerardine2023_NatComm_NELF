## This script establishes the environment used in the figure generation 
## notebooks. It relies on the data and other scripts present in this project
## folder.
##
## Running this entire script will load upwards of 3Gb of data into memory
## 

## ===OPTION 1===
## Set the number of CPU cores to be used for all computational steps
options(mc.cores = 8L)

## ==OPTION 2===
## Set the desired output folder for all generated plots. 
outdir <- here("figures/")


# Load Dependencies -------------------------------------------------------

suppressPackageStartupMessages({
    library(here)
    library(dplyr)
    library(reshape2)
    library(ComplexHeatmap)
    library(ggplot2)
    library(lemon) # vaporware ggplot2 extension; used for some facetting plots
    library(GGally)
    library(ggpubr)
    library(GenomicFeatures)
    library(DESeq2)
    library(BRGenomics)
    library(cetcolor)
    library(viridis)
    library(patchwork)
    library(BSgenome.Dmelanogaster.UCSC.dm6)
    require(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
    require(org.Dm.eg.db)
    library(Biostrings)
    library(ggseqlogo)
})

# Check user option
if (getOption("mc.cores") > detectCores()) {
    warning("Setting option \"mc.cores\" to not exceed number of virtual cores
            detected on your machine")
    options(mc.cores = detectCores())
}


# Source Scripts ----------------------------------------------------------

source(here("scripts/dataImport.R"))
source(here("scripts/utils.R"))
source(here("scripts/plot_shot_functions.R"))
source(here("scripts/plot_shot_wrapper.R"))


# Import Drosophila Annotations -------------------------------------------

dataImport.txs_all_usedregions()


# Import S2 Cell PRO-seq --------------------------------------------------

dataImport.PROseq.b1()
dataImport.PROseq.b2()


# Import PRO-seq read count tables (with spike-in NFs) --------------------

dataImport.readcounts_PROseq.b1()
dataImport.readcounts_PROseq.b2()



# Obtain compound NFs for PRO-seq (affecting FP conditions) ---------------

# Using long-gene-end normalization for FP conditions (and spike-in for DMSO),
# which only affects experiment/batch 1 (see the manuscript for info)

# I manually curated genes to drop from the long-gene end normalization.
#   These are based on regions in `txs.f_match` longer than 30 kb 
#   (so this is comprehensive to all regions that length or longer).
#   Genes are excluded for having obvious problems, like internal PRO-cap sites

names_drop <- c("zfh2", "Pde11", "ec", "sli", "luna", "B4", "CG32369", 
                "CG32982", "CG42684", "Ten-m", "mtgo", "robo2", "Lasp", 
                "pigs", "CG42674", "EcR", "InR", "CG34347", "Samuel", "Dyrk2", 
                "foxo", "tai", "skd", "N", "rgn", "ome", "NA", "LpR2", "Msp300",
                "plum", "scrib", "sdt", "MYPT-75D", "ttv", "dnc", "nmo", 
                "CG43367", "CG32521", "sdk", "htt", "Src64B", "norpA", "Pura",
                "Camta", "sm", "NetB", "CG45263", "Pvf3", "Hs6st", "CG13229", 
                "CG12541", "NetA", "axo", "CG5888", "Mp", "CG3655", 
                "Trim9", "CG42458", "ab", "Ac13E", "CG42673", "CG9650", "rg",
                "hiw", "CG43861", "beat-IIIc", "CG17716", "CG14431", "Myo10A",
                "stan", "DAAM", "grn")


# Get compound NFs (spike-in for DMSO, LGE for FP) for batch 1
subset(txs.f_match, width >= 5.5e4) %>% # more length-restricted set of genes
    subset(!symbol %in% names_drop) %>% 
    genebodies(5e4, -5e3) %>% 
    getCountsByRegions(PROseq.b1, .) %>% # get long gene-end counts
    lapply(sum) %>%
    as.data.frame %>% 
    melt(
        variable.name = "sample", 
        value.name = "lge_counts", 
        id.vars = NULL
    ) %>% 
    .$lge_counts %>% 
    
    (function(lge_counts) {
        
        # will modify a copy dataframe (little inelegant, copy-pasting...)
        df <- readcounts <- readcounts.b1
        df$KD <- sub("_(DMSO|FP)_rep.", "", df$sample)
        df$rep <- sub(".*_rep", "", df$sample)
        df$cond <- sub(".*(DMSO|FP).*", "\\1", df$sample)
        
        # get normalized long-gene-end counts (only need for DMSO conditions)
        norm_lgecounts <- lge_counts * df$NF
        
        # for FP, get NFs that match normalized long-gene-end counts in control
        readcounts$NF <- sapply(seq_len(nrow(df)), function(i) {
            # keep NFs for DMSO
            if (df$cond[i] == "DMSO") {
                return(df$NF[i])
            }
            
            # get within-rep/within-KD DMSO control
            idx_ctrl <- which(
                df$KD == df$KD[i] & 
                df$rep == df$rep[i] & 
                df$cond == "DMSO"
            )
            
            # get NF to equalize normalized long-gene-end readcounts
            round(norm_lgecounts[idx_ctrl] / lge_counts[i], 4)
        })
        
        readcounts
    }) -> readcounts.b1_alt



# Get compound NFs (spike-in for DMSO, LGE for FP) for batch 2
subset(txs.f_match, width >= 5.5e4) %>% # more length-restricted set of genes
    subset(!symbol %in% names_drop) %>% 
    genebodies(5e4, -5e3) %>% 
    getCountsByRegions(PROseq.b2, .) %>% # get long gene-end counts
    lapply(sum) %>%
    as.data.frame %>% 
    melt(
        variable.name = "sample", 
        value.name = "lge_counts", 
         id.vars = NULL
    ) %>% 
    
    .$lge_counts %>% 
    
    (function(lge_counts) {
        # more inelegant stuff from copy-pasting
        df <- readcounts <- readcounts.b2
        df$KD <- ifelse(grepl("LacZ", df$sample), "LacZ", "NELF")
        df$rep <- sub(".*_rep", "", df$sample)
        df$cond <- sub(".*(DMSO|FP_..min).*", "\\1", df$sample)
        
        # get normalized long-gene-end counts (only need for DMSO conditions)
        norm_lgecounts <- lge_counts * df$NF
        
        # for FP, get NFs that match normalized long-gene-end counts in control
        readcounts$NF <- sapply(seq_len(nrow(df)), function(i) {
            # keep NFs for DMSO
            if (df$cond[i] == "DMSO") {
                return(df$NF[i])
            }
            
            # get within-rep/within-KD DMSO control
            idx_ctrl <- which(
                df$KD == df$KD[i] & 
                df$rep == df$rep[i] & 
                df$cond == "DMSO"
            )
            
            # get NF to equalize normalized long-gene-end readcounts
            round(norm_lgecounts[idx_ctrl] / lge_counts[i], 4)
        })
        
        readcounts
    }) -> readcounts.b2_alt


# Import S2 Cell RNA-seq & NFs --------------------------------------------

dataImport.RNAseq()
dataImport.readcounts_RNAseq()



# Import imaging data -----------------------------------------------------

dataImport.PAGFP_imaging()
dataImport.MCPGFP_imaging()
dataImport.MCPGFP_imaging_prop_active()


invisible(gc())
