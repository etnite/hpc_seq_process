#!/bin/bash
shopt -s nullglob

## Interleave Paired Fastq Files
##
## This script is designed to take a single pair of fastq files, and interleave
## the reads into a single file using BBMap's reformat.sh.
##
## It takes a single string as its only positional argument. 
## Then it searches for corresponding fastq files using the pattern: 
## 
##     in_dir/<pattern>*R[1/2]*[fastq/fq].gz
##
## A search pattern can contain a subdirectory prefix. In this case, the
## subdirectory will be created in the output directory. 
##
## The output fastq file will be gzipped, consisting of the interleaved reads
## of the two input files. It will be written to the location:
##
##     out_dir/<pattern>_interleaved.fastq.gz
################################################################################


#### User-defined variables ####

## Input fastq directory, output directory for, and file
## of filename patterns, one listed per line
in_dir="/project/genolabswheatphg/merged_fastqs/SRW_excap"
out_dir="/project/genolabswheatphg/merged_fastqs/SRW_excap/fq_demux_test"
patterns_file="/home/brian.ward/test_samplist2.txt"


#### Executable ####

module load bbtools
module load pigz

echo
echo "interleave_fastqs.sh"
echo "Start time:"
date

mkdir -p "$out_dir"
array_ind=$1

## Get search pattern string and sample name; convert sample name to uppercase
patt=$(head -n "$array_ind" "$patterns_file" | tail -n 1)
sub_dir=$(dirname "$patt")
if [[ "$subdir" != "." ]]; then
	mkdir -p "${out_dir}/${sub_dir}"
fi

## Get the two paired fastq files
fq=$(echo "${in_dir}/${patt}"*R1*fastq.gz "${in_dir}/${patt}"*R1*fq.gz)
fq2=$(echo "${in_dir}/${patt}"*R2*fastq.gz "${in_dir}/${patt}"*R1*fq.gz)

## Interleave the files
reformat.sh in1="$fq" in2="$fq2" out="${out_dir}/${patt}_interleaved.fastq.gz"

echo
echo "End time:"
date
