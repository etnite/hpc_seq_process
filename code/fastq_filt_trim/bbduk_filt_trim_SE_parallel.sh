#!/bin/bash

## Single-end FASTQ Quality filtering and adapter trimming using BBDuk
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This version of the script runs BBDuk on a single fastq file.
## It takes a single string (usually a sample name) as its only positional
## argument. Then it searches for corresponding fastq files using the pattern:
## input_dir/samplename.fastq.gz
##
## This script is intended to be used with a parallel dispatch script, to enable independent
## parallel runs on multiple samples simultaneously.
##
## NOTES: 
##   1) bbduk parameters are hard coded in the bbduk.sh call below.
##   2) Working directory inherited from parallelizing script - it is easiest
##      to define absolute paths
################################################################################


#### User-defined constants ####

in_dir="/project/guedira_seq_map/Allegro_test/raw_fastq"
adapt_fasta="/home/brian.ward/Downloads/adapters.fa"
samp_file="/home/brian.ward/samp_and_file_lists/Allegro_raw_fq_samp_list.txt"
out_dir="/project/guedira_seq_map/Allegro_test/filt_fastq"


#### Executable  ####

module load singularity/3.7.1
module load bbmap

echo
echo "Start bbduk_filt_trim_paired_parallel.sh"
echo "Start time:"
date

mkdir -p "${out_dir}"
array_ind=$1

## Get length of adapters
#ad_len=$(head -n 2 "${adapt_fasta}" | tail -n -1 | wc -c)

## Get sample name
samp=$(head -n "${array_ind}" "${samp_file}" | tail -n 1)

## Run BBDuk
## Setting maq will remove entire reads based on their average qual
## maq=10 average error is 10%
## maq=13 ~5%
## maq=15 ~3%
## maq=20 1%
bbduk.sh -Xmx28g in="${in_dir}/${samp}.fastq.gz" \
    out="${out_dir}/${samp}.fastq.gz" \
    ref="${adapt_fasta}" ktrim=r k=25 mink=10 hdist=3 hdist2=1 ftm=5 maq=13 minlen=75

echo
echo "End time:"
date
