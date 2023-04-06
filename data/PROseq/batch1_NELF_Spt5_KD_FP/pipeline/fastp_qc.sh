#!/bin/bash

## This script generates fastp reports for each of the samples
## The motivation is for QC, and especially for checking insert size distributions
## This script is for single-ended sequencing

## This was made in hindsight, as the original processing pipeline did not retrieve
## insert sizes

## parallelization
let "ncores = $nproc - 2";

### This script is for UMIs (!) ###
# First 6 bases of read are UMI

RNA3p="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC";
nproc=`nproc`;

## script will make directories if they don't exist
[ -e logs/fastp_qc/insert_sizes ] || mkdir -p logs/fastp_qc/insert_sizes;
[ -e logs/fastp_qc/json ] || mkdir logs/fastp_qc/json;
[ -e logs/fastp_qc/html ] || mkdir logs/fastp_qc/html;

for fastq in "$@"
do

  outbase=$(echo $fastq | sed 's/.*\///' | sed 's/[.].*$//');

  [ -e tmp.fastq ] && rm tmp.fastq;

  echo " "
  echo "****************************************"
  echo "Begin: $outbase"
  echo "****************************************"
  echo " "

  fastp \
  --in1 $fastq \
  --out1 tmp.fastq \
  --adapter_sequence $RNA3p \
  --umi \
  --umi_loc=read1 \
  --umi_len=6 \
  --html logs/fastp_qc/html/${outbase}.html \
  --json logs/fastp_qc/json/${outbase}.json \
  --thread $ncores \
  --disable_length_filtering;

  cat tmp.fastq | awk '{if(NR%4==2) print length($1)}' | sort -n | uniq -c > logs/fastp_qc/insert_sizes/${outbase}.txt;
  rm tmp.fastq;
done;

exit 0;
