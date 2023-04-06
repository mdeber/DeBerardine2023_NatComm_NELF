#!/usr/bin/env bash

./pipeline/processPROseq \
--input-dir=bam \
--outdir-raw=bw/raw \
--outdir-norm=bw/norm \
--outdir-rds=rds_files \
--metadata-output=logs/readcounts.txt \
--spike-pattern=mm10 \
--control-pattern=LacZKD_DMSO \
--batch-norm=TRUE \
--paired-end=FALSE \
--quality=30 \
--threads=12 \
--threads-bam-import=10;
# --yield-size is NA (no chunking)
