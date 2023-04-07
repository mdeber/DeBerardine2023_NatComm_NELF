[![DOI](https://zenodo.org/badge/624615563.svg)](https://zenodo.org/badge/latestdoi/624615563)

# Overview

This repo contains all the code necessary to generate all figures in the
manuscript entitled "The NELF pausing checkpoint mediates the functional 
divergence of Cdk9" by Michael DeBerardine, Gregory T. Booth, Philip P. 
Versluis, and John T. Lis (Nature Communications 2023).

This repo also contains various ancillary data related to the sequencing 
experiments, as well as the tabulated live-imaging data used for plotting.

However, the sequencing data itself must be downloaded from GEO.

# Notes on sequencing data

Sequencing data must be downloaded from 
[GSE211397](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211397). 

All analysis code in this project utilizes the raw (not normalized) data (such
that all read counts are the original integer values), and normalization is
applied during each quantitative step.

However, we also recommend downloading the normalized data for the purposes of
genome browsing, e.g. using IGV.

# Instructions for downloading sequencing data

**STEP 1**

On the GEO page for 
[GSE211397](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211397), go to 
supplementary files and either download the entire archive, or using "custom", 
make sure all of the "raw" (non-normalized) files are selected. 

(You can optionally download the PROcap data files as well, but they aren't 
essential for generating the figure plots).


**STEP 2**

Extract the tar archive as a directory, which will be named `GSE211397_RAW`,
and put that directory into this parent directory 
(`DeBerardine2023_NatComm_NELF`).

**STEP 3**

You can now run the code in the `figure_generation.Rmd` script, which will
run various other scripts and functions to establish the interactive analysis 
environment in R and perform the various quantification and plotting 
operations.

# Other notes

The `data` directory did originally contain the fastq, bam, and bigWig/bedGraph 
files for each assay, but for simplicity, this repo keeps the relevant data 
files exactly as they're downloaded from GEO. However, the `data` directories 
still contain the pipeline files and logs, including the tabulations of the 
experimental vs. spike-in reads as well as the spike-in normalization factors. 

Note that for the Flavopiridol-treated PRO-seq samples, the "compound 
normalization factors" are generated in the `fig_generation_environment.R` 
script, while the normalization factors in the `readcounts.txt` files 
(found in subdirectories of the `data` directory) are the original (spike-in 
only) normalization factors.

# For Figure 5 (previous S. pombe data)

To run the code used to generate Figure 5, `GSE102308_COMBINED_bigwigs.tar.gz`
must be [downloaded from GEO](
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102308), extracted, and
the resulting directory must be placed in the `extdata` folder.

In `figure_generation.Rmd`, the code chunks for this figure are not run by 
default, and the data is not imported by default.

