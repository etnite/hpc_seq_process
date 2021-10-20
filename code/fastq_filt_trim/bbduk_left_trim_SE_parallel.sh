#!/bin/bash

## Single-end FASTQ lefthand trim using BBDuk
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script performs a simple constant left-hand trim of all reads within a
## .fastq file. This is useful when we expect the beginning of each read to
## contain, for instance, a primer or probe which may not contain polymorphisms
## that are actually present in its binding site.
##
## This script is intended to be used with a parallel dispatch script, to enable independent
## parallel runs on multiple samples simultaneously.
##
## The user needs to define the number of bases to trim from the left-hand side
## of each read (trim_size) and the output directory. The list of input files
## and the index specifying which file from the list to process are inherited
## from the parallel dispatch script.
##
## NOTES:
##   1) bbduk parameters are hard coded in the bbduk.sh call below.
##   2) Working directory inherited from parallelizing script - it is easiest
##      to define absolute paths
################################################################################


#### User-defined constants ####

trim_size=40
out_dir="/project/guedira_seq_map/Allegro_test/left_trim_fastq"


#### Executable  ####

module load singularity/3.7.1
module load bbmap

echo
echo "Start bbduk_filt_trim_paired_parallel.sh"
echo "Start time:"
date

mkdir -p "${out_dir}"
fasta_list=$1
array_ind=$2

## Get length of adapters
#ad_len=$(head -n 2 "${adapt_fasta}" | tail -n -1 | wc -c)

## Get sample name
infile=$(head -n "${array_ind}" "${fasta_list}" | tail -n 1)
base=$(basename "$infile")
outfile="${out_dir}/${base}"

## Run BBDuk
bbduk.sh -Xmx28g in="$infile" out="${out_dir}/${base}" ftl=$trim_size

echo
echo "End time:"
date
