#!/usr/bin/bash

# calling STAR 2.7.3a using ENCODE settings, as found in the documentation

STAR \
--genomeDir dm6_mm10_stargenome \
--genomeLoad LoadAndExit;

for file in "$@";
do

outbase=$(echo $file | sed 's:.*/::' | sed 's/.fastq.gz//' | sed 's/trim_//');

[ -e logs/umi_tools ] || mkdir logs/umi_tools

# have given STAR limit of 12Gb for sorting bam file
STAR \
--runThreadN 12 \
--genomeDir dm6_mm10_stargenome \
--genomeLoad LoadAndKeep \
--readFilesIn $file \
--readFilesCommand zcat \
--outFileNamePrefix bam/${outbase} \
--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--alignIntronMin 20 \
--alignIntronMax 1000000;

# rename the idiotic and stupid output names STAR forces on you
mv bam/${outbase}Aligned.out.sam bam/${outbase}.sam

# sort and index bam and remove PCR duplicates
samtools view -b -S bam/${outbase}.sam |
samtools sort -@ 12 -o bam/${outbase}.bam;
samtools index bam/${outbase}.bam;

# fastp using ":" to separater, umi_tools default is "_"
umi_tools dedup \
-I bam/${outbase}.bam \
--umi-separator=: \
--output-stats=logs/umi_tools/dedupstats_${outbase}.tsv \
-S bam/${outbase}.bam;

rm bam/{outbase}.sam;

done;

STAR \
--genomeDir dm6_mm10_stargenome \
--genomeLoad Remove;

exit 0;
