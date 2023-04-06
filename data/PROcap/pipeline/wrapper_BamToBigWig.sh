#!/bin/bash
Rscript pipeline/BamToBigWig_v2.R \
-i bam \
-o bw \
--genome dm6 \
--revcomp F \
-e 5p \
--stranded T \
--spike-genome mm10 \
--quality 30 \
--threads 12;
