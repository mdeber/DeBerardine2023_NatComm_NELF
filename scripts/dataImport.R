
suppressPackageStartupMessages({
    require(here)
    require(GenomicFeatures)
    require(BRGenomics)
    require(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
    require(org.Dm.eg.db)
})


# Data import helper function ---------------------------------------------

.dataImport.checkFiles <- function(assay = NULL, batch = NULL) {
    # Called by the data import function, this helper function checks that the 
    # required files exist
    
    if (!"GSE211397_RAW" %in% dir(here())) {
        stop("Please see the instructions in the README.md file. Data from the
        GSE211397_RAW tarball must be downloaded from GEO, extracted
        as a directory, and placed within DeBerardine2023_NatComm_NELF before
        any data can be imported."
        )
    }
    
    files_proseq_b1 <- c(
        "GSM6469458_PROseq_exp1_LacZKD_DMSO_rep1_raw_minus.bw",
        "GSM6469458_PROseq_exp1_LacZKD_DMSO_rep1_raw_plus.bw",
        "GSM6469459_PROseq_exp1_LacZKD_DMSO_rep2_raw_minus.bw",
        "GSM6469459_PROseq_exp1_LacZKD_DMSO_rep2_raw_plus.bw",
        "GSM6469460_PROseq_exp1_LacZKD_FP_rep1_raw_minus.bw",
        "GSM6469460_PROseq_exp1_LacZKD_FP_rep1_raw_plus.bw",
        "GSM6469461_PROseq_exp1_LacZKD_FP_rep2_raw_minus.bw",
        "GSM6469461_PROseq_exp1_LacZKD_FP_rep2_raw_plus.bw",
        "GSM6469462_PROseq_exp1_NELFeKD_DMSO_rep1_raw_minus.bw",
        "GSM6469462_PROseq_exp1_NELFeKD_DMSO_rep1_raw_plus.bw",
        "GSM6469463_PROseq_exp1_NELFeKD_DMSO_rep2_raw_minus.bw",
        "GSM6469463_PROseq_exp1_NELFeKD_DMSO_rep2_raw_plus.bw",
        "GSM6469464_PROseq_exp1_NELFeKD_FP_rep1_raw_minus.bw",
        "GSM6469464_PROseq_exp1_NELFeKD_FP_rep1_raw_plus.bw",
        "GSM6469465_PROseq_exp1_NELFeKD_FP_rep2_raw_minus.bw",
        "GSM6469465_PROseq_exp1_NELFeKD_FP_rep2_raw_plus.bw"
    )
    
    files_proseq_b2 <- c(
        "GSM6469466_PROseq_exp2_LacZKD_DMSO_rep1_raw_minus.bw",
        "GSM6469466_PROseq_exp2_LacZKD_DMSO_rep1_raw_plus.bw",
        "GSM6469467_PROseq_exp2_LacZKD_DMSO_rep2_raw_minus.bw",
        "GSM6469467_PROseq_exp2_LacZKD_DMSO_rep2_raw_plus.bw",
        "GSM6469468_PROseq_exp2_LacZKD_FP_05min_rep1_raw_minus.bw",
        "GSM6469468_PROseq_exp2_LacZKD_FP_05min_rep1_raw_plus.bw",
        "GSM6469469_PROseq_exp2_LacZKD_FP_05min_rep2_raw_minus.bw",
        "GSM6469469_PROseq_exp2_LacZKD_FP_05min_rep2_raw_plus.bw",
        "GSM6469470_PROseq_exp2_LacZKD_FP_10min_rep1_raw_minus.bw",
        "GSM6469470_PROseq_exp2_LacZKD_FP_10min_rep1_raw_plus.bw",
        "GSM6469471_PROseq_exp2_LacZKD_FP_10min_rep2_raw_minus.bw",
        "GSM6469471_PROseq_exp2_LacZKD_FP_10min_rep2_raw_plus.bw",
        "GSM6469472_PROseq_exp2_LacZKD_FP_20min_rep1_raw_minus.bw",
        "GSM6469472_PROseq_exp2_LacZKD_FP_20min_rep1_raw_plus.bw",
        "GSM6469473_PROseq_exp2_LacZKD_FP_20min_rep2_raw_minus.bw",
        "GSM6469473_PROseq_exp2_LacZKD_FP_20min_rep2_raw_plus.bw",
        "GSM6469474_PROseq_exp2_NELFeKD_DMSO_rep1_raw_minus.bw",
        "GSM6469474_PROseq_exp2_NELFeKD_DMSO_rep1_raw_plus.bw",
        "GSM6469475_PROseq_exp2_NELFeKD_DMSO_rep2_raw_minus.bw",
        "GSM6469475_PROseq_exp2_NELFeKD_DMSO_rep2_raw_plus.bw",
        "GSM6469476_PROseq_exp2_NELFeKD_FP_05min_rep1_raw_minus.bw",
        "GSM6469476_PROseq_exp2_NELFeKD_FP_05min_rep1_raw_plus.bw",
        "GSM6469477_PROseq_exp2_NELFeKD_FP_05min_rep2_raw_minus.bw",
        "GSM6469477_PROseq_exp2_NELFeKD_FP_05min_rep2_raw_plus.bw",
        "GSM6469478_PROseq_exp2_NELFeKD_FP_10min_rep1_raw_minus.bw",
        "GSM6469478_PROseq_exp2_NELFeKD_FP_10min_rep1_raw_plus.bw",
        "GSM6469479_PROseq_exp2_NELFeKD_FP_10min_rep2_raw_minus.bw",
        "GSM6469479_PROseq_exp2_NELFeKD_FP_10min_rep2_raw_plus.bw",
        "GSM6469480_PROseq_exp2_NELFeKD_FP_20min_rep1_raw_minus.bw",
        "GSM6469480_PROseq_exp2_NELFeKD_FP_20min_rep1_raw_plus.bw",
        "GSM6469481_PROseq_exp2_NELFeKD_FP_20min_rep2_raw_minus.bw",
        "GSM6469481_PROseq_exp2_NELFeKD_FP_20min_rep2_raw_plus.bw"
    )
    
    files_rnaseq <- c(
        "GSM6469482_RNAseq_exp1_LacZKD_DMSO_rep1_raw_minus.bedGraph.gz",
        "GSM6469482_RNAseq_exp1_LacZKD_DMSO_rep1_raw_plus.bedGraph.gz",
        "GSM6469483_RNAseq_exp1_LacZKD_DMSO_rep2_raw_minus.bedGraph.gz",
        "GSM6469483_RNAseq_exp1_LacZKD_DMSO_rep2_raw_plus.bedGraph.gz",
        "GSM6469484_RNAseq_exp1_NELFeKD_DMSO_rep1_raw_minus.bedGraph.gz",
        "GSM6469484_RNAseq_exp1_NELFeKD_DMSO_rep1_raw_plus.bedGraph.gz",
        "GSM6469485_RNAseq_exp1_NELFeKD_DMSO_rep2_raw_minus.bedGraph.gz",
        "GSM6469485_RNAseq_exp1_NELFeKD_DMSO_rep2_raw_plus.bedGraph.gz"
    )
    
    files_procap <- c(
        "GSM6469486_PROcap_LacZKD_minus.bw",
        "GSM6469486_PROcap_LacZKD_plus.bw",
        "GSM6469487_PROcap_LacZKD_noTAP_minus.bw",
        "GSM6469487_PROcap_LacZKD_noTAP_plus.bw"
    )
    
    
    
    if (is.null(assay) & is.null(batch)) {
        req_files <- c(
            files_proseq_b1,
            files_proseq_b2,
            files_rnaseq,
            files_procap
        )
    } else if (assay == "PROseq") {
        if (is.null(batch)) {
            req_files <- c(
                files_proseq_b1,
                files_proseq_b2
            )
        } else if (batch == 1) {
            req_files <- files_proseq_b1
        } else if (batch == 2) {
            req_files <- files_proseq_b2
        } else {
            stop("for assay=PROseq, batch must be NULL, 1, or 2 (numeric)")
        }
    } else if (assay == "RNAseq") {
        req_files <- files_rnaseq
    } else if (assay == "PROcap") {
        req_files <- files_procap
    } else {
        stop("assay must be one of \"PROseq\", \"RNAseq\", or \"PROcap\"")
    }
    
    
    if (!all(req_files %in% dir(here("GSE211397_RAW")))) {
        stop("You've attempted to import data that isn't found within 
        DeBerardine2023_NatComm_NELF/GSE211397_RAW. For the given assay, make
        sure all non-normalized files from GSE211397 are in the directory
        DeBerardine2023_NatComm_NELF/GSE211397_RAW, and do not rename them. 
        See README.md for instructions")
    }
}


# PRO-seq import functions ------------------------------------------------

dataImport.PROseq.b1 <- function() {
    # Import batch 1 PROseq data as a list of GRanges objects and assign that in 
    # the global environment
    
    .dataImport.checkFiles(assay = "PROseq", batch = 1)
    
    t.start <- Sys.time()
    PROseq.b1 <<- here("GSE211397_RAW") %>% 
        dir(full.names = TRUE) %>%
        .[grep("PROseq_exp1(.*)_raw_", .)] %>% 
        split(., grepl("minus.bw", .)) %>% 
        setNames(c("plus_file", "minus_file")) %>% 
        c(genome = "dm6") %>% 
        do.call(import_bigWig, .) %>% 
        setNames(., sub(".*PROseq_exp1_(.*)_raw.*", "\\1", names(.)))
    
    t.diff <- Sys.time() - t.start
    cat(
        "Imported PROseq.b1 in", 
        round(t.diff[[1]], 1), 
        units(t.diff), 
        "\n"
    )
}

dataImport.PROseq.b2 <- function() {
    # Import batch 2 PROseq data as a list of GRanges objects and assign that in 
    # the global environment
    
    .dataImport.checkFiles(assay = "PROseq", batch = 2)
    
    t.start <- Sys.time()
    PROseq.b2 <<- here("GSE211397_RAW") %>% 
        dir(full.names = TRUE) %>%
        .[grep("PROseq_exp2(.*)_raw_", .)] %>% 
        split(., grepl("minus.bw", .)) %>% 
        setNames(c("plus_file", "minus_file")) %>% 
        c(genome = "dm6") %>% 
        do.call(import_bigWig, .) %>% 
        setNames(., sub(".*PROseq_exp2_(.*)_raw.*", "\\1", names(.)))
    t.diff <- Sys.time() - t.start
    cat(
        "Imported PROseq.b2 in", 
        round(t.diff[[1]], 1), 
        units(t.diff), 
        "\n"
    )
}



# RNA-seq import function -------------------------------------------------

dataImport.RNAseq <- function() {
    # Import RNAseq data (batch 1) as a list of GRanges objects and assign that 
    # in  the global environment
    .dataImport.checkFiles(assay = "RNAseq")
    
    t.start <- Sys.time()
    
    RNAseq.b1 <<- here("GSE211397_RAW") %>% 
        dir(full.names = TRUE) %>%
        .[grep("RNAseq_exp1(.*)_raw_", .)] %>% 
        split(., grepl("minus.bedGraph.gz", .)) %>% 
        setNames(c("plus_file", "minus_file")) %>% 
        c(genome = "dm6") %>% 
        do.call(import_bedGraph, .) %>% 
        setNames(., sub(".*RNAseq_exp1_(.*)_raw.*", "\\1", names(.)))
    
    t.diff <- Sys.time() - t.start
    cat(
        "Imported RNAseq.b1 in", 
        round(t.diff[[1]], 1), 
        units(t.diff), 
        "\n"
    )
}


# PRO-cap import function -------------------------------------------------

dataImport.PROcap <- function() {
    # Import PROcap data as a list of GRanges objects and assign that in the 
    # global environment
    .dataImport.checkFiles(assay = "PROcap")
    
    t.start <- Sys.time()
    
    PROcap <<- here("GSE211397_RAW") %>% 
        dir(full.names = TRUE) %>%
        .[grep("PROcap", .)] %>% 
        split(., grepl("minus.bw", .)) %>% 
        setNames(c("plus_file", "minus_file")) %>% 
        c(genome = "dm6") %>% 
        do.call(import_bigWig, .) %>% 
        setNames(., sub(".*PROcap_(.*)_(minus|plus).bw", "\\1", names(.)))
    
    t.diff <- Sys.time() - t.start
    cat(
        "Imported PROcap in", 
        round(t.diff[[1]], 1), 
        units(t.diff), 
        "\n"
    )
}


# Readcount table import functions ----------------------------------------

dataImport.readcounts_PROseq.b1 <- function() {
    # Import read counts table for batch 1 (experiment 1) PRO-seq and assign it
    # in the global environment
    readcounts.b1 <<- read.delim(
        here("data/PROseq/batch1_NELF_Spt5_KD_FP/logs/readcounts.txt"),
        stringsAsFactors = FALSE
    )
    cat("Imported readcounts.b1", "\n")
}

dataImport.readcounts_PROseq.b2 <- function() {
    # Import read counts table for batch 2 (experiment 2) PRO-seq and assign it
    # in the global environment
    readcounts.b2 <<- read.delim(
        here("data/PROseq/batch2_NELF_KD_FPTC/logs/readcounts.txt"),
        stringsAsFactors = FALSE
    )
    cat("Imported readcounts.b2", "\n")
}


dataImport.readcounts_RNAseq <- function() {
    # Import read counts table for RNA-seq (batch 1/experiment 1) and assign it
    # in the global environment
    readcounts_RNAseq.b1 <<- read.delim(
        here("data/RNAseq/logs/readcounts.txt"),
        stringsAsFactors = FALSE
    )
    cat("Imported readcounts_RNAseq.b1", "\n")
}


# Annotation Import Functions ---------------------------------------------

dataImport.txdb_txs <- function() {
    ## Makes a GRanges object with all annotated transcripts, and also uses
    ## as org.db object to get the fly gene names (symbols)
    txdb <- 
        TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene
    
    x <- tidyChromosomes(transcripts(txdb))
    
    suppressMessages({
        x$gene_id <- select(
            txdb, 
            keys = x$tx_name, 
            keytype = "TXNAME", 
            columns = "GENEID"
        )[[2]]
        
        x$symbol <-select(
            org.Dm.eg.db::org.Dm.eg.db, 
            keys = x$gene_id,
            keytype = "FLYBASE", 
            columns = "SYMBOL"
        )[[2]]
    })
    x$tx_id <- NULL
    x
}


dataImport.txs_5p <- function() {
    # PRO-cap filtered TSSes, but gene bodies are not filtered at all.
    #   (so could overlap other TSSes, multiple CPSes, etc.)
    #   -> so this is only used for its promoter regions
    # Briefly, takes most 5' isoform and keeps only if >40% of PRO-cap in 
    #   [50, +150] is found in that TSS.
    
    txs.5p <- import.bed(here("annotations/filtered_5p_only_txs.bed"))
    mcols(txs.5p) <- read.delim(
        here("annotations/filtered_5p_only_txs_names.tsv"),
        stringsAsFactors = FALSE
    )
    txs.5p
}


dataImport.txs_filt <- function() {
    ## Heavily-filtered transcripts, with isoform selection, earliest CPS 
    ## selection, PRO-cap concentrated on a single site, no upstream readthrough
    txs <- import.bed(here("annotations/filtered_txs.bed"))
    names(mcols(txs)) <- c("gene_id", "PROcap")
    
    txs_names <- read.delim(
        here("annotations/filtered_txs_names.tsv"),
        stringsAsFactors = FALSE
    )
    
    # will throw an error if they're not the same
    txs$symbol <- txs_names$symbol[txs_names$gene_id == txs$gene_id]
    txs
}


dataImport.txs_gb_consensus <- function() {
    ## conservative ("central" or "consensus") gene bodies
    ##   most downstream TSS + 300 bp
    ##   most upstream CPS - 300 bp
    ##   of that, removed overlapping gene-bodies
    ## not suitable for anything to do with pausing, but adds a lot more genes 
    ## for evaluating expression (gene body Pol II density)
    
    txs <- import.bed(here("annotations/consensus_gb.bed"))
    names(mcols(txs)) <- c("gene_id", "PROseq")
    
    txs_names <- read.delim(
        here("annotations/consensus_gb_names.tsv"),
        stringsAsFactors = FALSE
    )
    
    # will throw an error if they're not the same
    txs$symbol <- txs_names$symbol[txs_names$gene_id == txs$gene_id]
    txs
}


dataImport.txs_cps <- function() {
    txs.cps <- import.bed(here("annotations/annotated_txends.bed"))
    mcols(txs.cps) <- read.delim(
        here("annotations/annotated_txends_names.tsv"),
        stringsAsFactors = FALSE
    )
    txs.cps
}


dataImport.txs_all_usedregions <- function() {
    ## this function does the main extra filtering steps that I do 
    ## to obtain the various final gene lists
    
    # first, import a txdb object, and all annotated transcripts
    txdb <<- 
        TxDb.Dmelanogaster.UCSC.dm6.ensGene::TxDb.Dmelanogaster.UCSC.dm6.ensGene
    cat("Imported txdb (dm6.ensGene)", "\n") 
    txs <<- dataImport.txdb_txs()
    cat("Imported txs (all annotated dm6 transcripts)", "\n") 
    
    # then my filtered gene lists
    txs.gb <<- dataImport.txs_gb_consensus()
    cat("Imported txs.gb", "\n") 
    
    # txs.5p = filtered promoter regions (from 5' transcripts)
    #   but gene bodies unsafe (not filtered)
    txs.5p <<- dataImport.txs_5p()
    txs.5p_pp <<- promoters(txs.5p, 0, 100)
    cat("Imported txs.5p", "\n") 
    cat("Imported txs.5p_pp", "\n")
    
    txs.f <<- dataImport.txs_filt()
    cat("Imported txs.f", "\n") 
    
    # get gene lists further filtered for matched promoters & gene bodies
    txs.f_pp_all <<- promoters(txs.f, 0, 100)
    txs.f_gb <<- genebodies(txs.f, 300, -300, min.window = 500)
    cat("Imported txs.f_pp_all", "\n")
    cat("Imported txs.f_gb", "\n")
    
    # promoter regions matching genes in txs.f_gb
    txs.f_pp_match_gb <<- subset(txs.f_pp_all, gene_id %in% txs.f_gb$gene_id)
    cat("Imported txs.f_pp_match_gb", "\n")
    
    # long genes, matching promoters & DISTAL gene bodies; (only genes > 6.2kb)
    txs.f_gb.long <<- genebodies(txs.f_gb, 5700, 0, min.window = 500)
    txs.f_pp_match_gb.long <<- subset(
        txs.f_pp_match_gb, 
        gene_id %in% txs.f_gb.long$gene_id
    )
    
    cat("Imported txs.f_gb.long", "\n")
    cat("Imported txs.f_pp_match_gb.long", "\n")
    
    # full genes for genes in matched promoter/gene body lists
    txs.f_match <<- subset(txs.f, gene_id %in% txs.f_gb$gene_id)
    txs.f_match.long <<- subset(txs.f, gene_id %in% txs.f_gb.long$gene_id)
    
    cat("Imported txs.f_match", "\n")
    cat("Imported txs.f_match.long", "\n")
    
    # txs.cps regions, but only those matching txs.gb genes;
    #   and also match those regions in txs.gb
    txs.cps <- dataImport.txs_cps()
    txs.cps <<- subset(txs.cps, gene_id %in% txs.gb$gene_id)
    txs.gb_match_cps <<- subset(txs.gb, gene_id %in% txs.cps$gene_id)
    
    cat("Imported txs.cps", "\n")
    cat("Imported txs.gb_match_cps", "\n")
}




# Import S. pombe PRO-seq data --------------------------------------------


dataImport.spomProseq <- function() {
    # Importing the S. pombe Cdk9as inhibition time-course PRO-seq data from
    # Booth et al. 2018. Calling it "batch 2" because there was another set of
    # experiments on S. pombe Cdk9as before the time-course.
    # 
    # The data is already spike-in normalized, although not in "portable" units
    # (like the RRPM I use for my data)
    
    t.start <- Sys.time()
    
    spomPROseq.b2_combined <<- here("extdata/GSE102308_COMBINED_bigwigs") %>% 
        dir(full.names = TRUE) %>% 
        .[grepl("ePROseq", .) & !grepl("Lsk1as", .)] %>% 
        split(., grepl("minus.bw", .)) %>% 
        setNames(c("plus_file", "minus_file")) %>% 
        do.call(import_bigWig, .) %>% 
        setNames(., {
            names(.) %>% 
                sub(".*ePROseq_Pombe_", "", .) %>% 
                sub("_COMBINED_pombe.*", "", .) %>% 
                sub("_(3MBPP1_10uM|DMSO)", "", .) %>% 
                sub("_(.)min", "_0\\1min", .) %>% 
                sub("_30sec", "_00min30sec", .)
        }) %>% 
        .[order(names(.))] %>% 
        
        # bigWig files have "standard" S. pombe chromosome names (Roman 
        # numerals), but the annotations I have from the Booth et al. 2018
        # manuscript use UCSC-style numbers. Will change the data to match the
        # annotations
        (function(x) {
            # will write them plainly, so first check that they are as
            # expected
            stopifnot(
                lapply(x, seqlevels) %>% 
                    unname %>% 
                    (function(x) sapply(x[-1], "==", x[[1]])) %>% 
                    all
            )
            x
        }) %>% 
        lapply(function(x) {
            # write the new chromosome names
            seqlevels(x) <- c(
                "chrAB325691", "chr01", "chr02", "chr03", "chrMT", "chrMTR"
            )
            x
        })
    
    t.diff <- Sys.time() - t.start
    cat(
        "Imported spomPROseq.b2_combined in", 
        round(t.diff[[1]], 1), 
        units(t.diff), 
        "\n"
    )
}



# Import S. pombe annotations ---------------------------------------------

dataImport.spom_txs_filt <- function() {
    # The filtered gene list from the Booth et al. 2018 paper ()
    spom_txs_filt <<- import.bed(
        here("extdata/PombeInfo_Booth2018/SP_CompleteFilteredGenes.bed")
    )
    
    cat("Imported spom_txs_filt", "\n")
    
    # Check that we get the 42 filtered genes >6kb expected from the paper
    stopifnot(
        sum(width(spom_txs_filt) > 6e3) == 42
    )
}



# Import Drosophila salivary gland imaging data ---------------------------

dataImport.PAGFP_imaging <- function() {
    # Import Rpb9::Photoactivatable_GFP imaging data tables
    
    Imaging_Rpb9PA <<- 
        here(
            "data/Imaging", 
            "Rpb9-PA-GFP_datatables"
        ) %>% 
        dir(full.names = TRUE) %>% 
        (function(x) {
            setNames(
                lapply(x, read.delim),
                paste0("x", sub(".*/", "", x))
            )
        }) %>% 
        dfList2df %>% 
        mutate(
            sample = sub("^x", "", sample),
            sample = sub(".txt", "", sample),
            kd = sub(".*(Lucif|NELF).*", "\\1", sample),
            kd = sub("Lucif", "Luciferase", kd),
            kd = paste(kd, "RNAi"),
            drug = sub(".*_(noFP|FP)_.*", "\\1", sample),
            drug = ifelse(drug == "noFP", "-FP", "+FP"),
            condition = paste(kd, drug),
            date = sub("_.*", "", sample)
        )
    
    cat("Imported Imaging_Rpb9PA", "\n")
}


dataImport.MCPGFP_imaging <- function() {
    # Import imaging data for MCP::GFP (intensity vs. background)
    
    Imaging_MCPGFP_intensity <<- 
        here(
            "data/Imaging", 
            "MCP-GFP_summary_datatables",
            "intensity_vs_background"
        ) %>% 
        dir(full.names = TRUE) %>% 
        (function(x) {
            setNames(
                lapply(x, read.delim, header = FALSE),
                sub(".*/(.*).txt", "\\1", x)
            )
        }) %>% 
        lapply(setNames, "MCPGFP_vs_background") %>% 
        dfList2df %>% 
        mutate(
            kd = sub(".*(Lucif|NELF).*", "\\1", sample),
            kd = sub("Lucif", "Luciferase", kd),
            kd = paste(kd, "RNAi"),
            drug = sub(".*_(noFP|FP)", "\\1", sample),
            drug = ifelse(drug == "noFP", "-FP", "+FP"),
            condition = paste(kd, drug)
        )
    
    cat("Imported Imaging_MCPGFP_intensity", "\n")
}


dataImport.MCPGFP_imaging_prop_active <- function() {
    # Import tabulation of active loci according to MCP::GFP imaging
    
    Imaging_MCPGFP_prop_active <<- 
        here(
            "data/Imaging",
            "MCP-GFP_summary_datatables",
            "proportion_loci_active.txt"
        ) %>% 
        read.delim(header = FALSE, sep = "\t") %>% 
        setNames(c("sample", "MCPGFP_prop_active")) %>% 
        mutate(
            kd = sub(".*(Lucif|NELF).*", "\\1", sample),
            kd = sub("Lucif", "Luciferase", kd),
            kd = paste(kd, "RNAi"),
            drug = sub(".*_(noFP|FP)", "\\1", sample),
            drug = ifelse(drug == "noFP", "-FP", "+FP"),
            condition = paste(kd, drug)
        )
    
    cat("Imported Imaging_MCPGFP_prop_active", "\n")
}





