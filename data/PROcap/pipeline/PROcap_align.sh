#!/bin/bash
## USAGE:  Provide a list of single-ended fastq/fastq.gz files to be aligned
## Example: bash ./PROseq_align.sh ./data/PROseq*R1.fastq.gz
## Script will trim 3' adapters, deplete rDNA, and aligns to combined genome

## Makes us of a masked combined genome, which lacks loci in mouse genome to
##  which Drosophila reads can map

####################################################################
### In this script, removed length trimming & length requirement ###
####################################################################

# script for UMIs - uses prinseq_lite
# prinseq_lite="/home/macproadmin/prinseq-lite-0.20.4/prinseq-lite.pl"

# cutadapt can only multithread if installed with python3.3+
nproc=`nproc`; # nproc command for Debian/Ubuntu

# path to bowtie2 genome indices
genome_mask="/home/macproadmin/users/MikeD/Mask_mm10/dm6_mm10_masked_bt2/dm6_mm10_masked";
rdna_dmel="/home/macproadmin/users/MikeD/rDNA_repeats/rDNAdmel_bt2/rDNAdmel";
rdna_mmus="/home/macproadmin/users/MikeD/rDNA_repeats/rDNAmouse_bt2/rDNAmouse";

# 3' adapter is NNNNNNGAUCGUCGGACUGUAGAACUCUGAAC
# so in read: GTTCAGAGTTCTACAGTCCGACGATCNNNNNN

# 3' adapter, for any remaining TRU-seq
RNA3p="TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC";

for fastq in "$@"
do
  if [ -z `basename $fastq | grep -i .fastq` ]
  then
      echo `basename $fastq` "does not have .fastq suffix - aborting";
      exit 1;
  fi
done

for fastq in "$@"
do
  dname=`dirname $fastq`; # output directory
  fname=`basename $fastq`; # output file
  fpath=$dname/${fname%.fastq*}; # file path
  gname=`basename $genome`; # genome name
  ext=${fname#*.}; # extension of input files
  # outbase=$dname/${gname}_spikecount/${fname%_R1.fastq*}; # base filename for output
  outbase=${fname%.fastq*};

  echo " "
  echo "****************************************"
  echo "Begin: $fname"
  echo "****************************************"
  echo " "

  # will use named pipes to avoid read/write bottlenecks
  # will check that they don't exist first
  [ -e /tmp/dep_dmel_rdna ] && rm /tmp/dep_dmel_rdna;
  [ -e /tmp/dep_mmus_rdna ] && rm /tmp/dep_mmus_rdna;
  [ -e /tmp/alignment_out ] && rm /tmp/alignment_out;
  mkfifo /tmp/dep_dmel_rdna;
  mkfifo /tmp/dep_mmus_rdna;
  mkfifo /tmp/alignment_out;

  # ADAPTER TRIMMING (cutadapt)
  ## --nextseq-trim is --quality-cutoff that also trims terminal high quality G's
  ## --cut 6 removes the first 6 bases
  cutadapt \
  -a $RNA3p \
  --cores $nproc \
  --nextseq-trim 10 \
  --minimum-length 15 \
  ${fpath}.$ext |

  # Align to Dmel rDNA repeat & keep unmapped
  bowtie2 \
  --very-sensitive \
  -x $rdna_dmel \
  -U /dev/stdin \
  -S /dev/null \
  --un /tmp/dep_dmel_rdna &

  # Align to Mmus rDNA repeat & keep unmapped
  bowtie2 \
  --very-sensitive \
  -x $rdna_mmus \
  -U /tmp/dep_dmel_rdna \
  -S /dev/null \
  --un /tmp/dep_mmus_rdna &

  # Align to combined genome with mm10 masked for Dmel PRO-seq alignments
  bowtie2 \
  --very-sensitive \
  --threads $nproc \
  -x $genome_mask \
  -U /tmp/dep_mmus_rdna \
  -S /tmp/alignment_out &

  # sort by name & convert to BAM
  samtools sort -@ $nproc -n -o ${outbase}.bam /tmp/alignment_out;

  rm /tmp/dep_dmel_rdna;
  rm /tmp/dep_mmus_rdna;
  rm /tmp/alignment_out;

  [ -e ${outbase}.bam ] || exit 1; # exit with error if file not made

  echo " "
  echo "Finished: $fname"
  echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  echo " "
done
