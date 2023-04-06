#!/bin/bash

# makes bedGraphs of all reads (identical reads are single lines with score = 2)

Rscript pipeline/process3pRNAseq.R \
--input-dir=bam_b1 \
--outdir-raw=bedgraph/raw \
--outdir-norm=bedgraph/norm \
--outdir-rds=rds_files_b1 \
--metadata-output=logs/readcounts_b1.tsv \
--spike-pattern=mm10 \
--control-pattern=b1_LacZ \
--batch-norm=TRUE \
--quality 30 \
--threads 12;
