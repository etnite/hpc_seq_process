#!/bin/bash

## Concatenate paired fastq files
##
## This script takes as input a single string (typically a genotype NAME). It
## will typically be used with arrayer.sh to run in independent, parallel
## fashion.
##
## The user must define an input directory, and an output directory.
##
## The fastq files should be in the input directory, named in the format 
## *NAME*_R1_*.fastq.gz for forward reads and *NAME*_R2_*.fastq.gz for reverse reads.
##
## The script will match any fastq files matching these patterns, and concatenate
## them into a single fastq file. It doesn't matter if the files are gzipped or
## not.
##
## It will also convert all sample names to uppercase, and convert underscores
## within sample names to dashes.
##
## NOTES:
##   1) THIS SCRIPT DOES NOT CREATE A SINGLE INTERLEAVED FASTQ - instead it
##      concatenates all forward reads into one file, and all reverse reads into
##      a second file
##   2) Working directory inherited from parallelizing script - it is easiest
##      to define absolute paths
################################################################################


#### User-defined variables ####

## No trailing slashes!
#in_dir="/project/genolabswheatphg/raw_data/Mary_untarred"
#out_dir="/project/genolabswheatphg/merged_fastqs/KS_HRW_excap"
#samp_file="/home/brian.ward/repos/wheat_phg/sample_lists/KS_HRW_samples.txt"

#in_dir="/project/genolabswheatphg/raw_data/wheatCAP_parents/separated_fq"
#out_dir="/project/genolabswheatphg/merged_fastqs/wheatCAP_parents"
#samp_file="/home/brian.ward/repos/wheat_phg/sample_lists/wheatCAP_samples.txt"

in_dir="/project/genolabswheatphg/raw_data/v1_hapmap"
out_dir="/project/genolabswheatphg/merged_fastqs/v1_hapmap"
samp_file="/home/brian.ward/repos/wheat_phg/sample_lists/v1_hapmap_bioproj/sample_names.txt"


#### Executable ####

echo
echo "Start concat_fastqs.sh"
echo "Start time:"
date

array_ind=$1
mkdir -p "${out_dir}"

## Get sample name
name=$(head -n "${array_ind}" "${samp_file}" | tail -n 1)

## Enable recursive globbing
shopt -s globstar nullglob

## Convert sample name to uppercase and underscores to dashes
## Personal pref - I like samples to have dashes to make pattern recognition easier
upname="${name^^}"
upname=$(echo "${upname}" | sed 's/_/-/g')

## Concatenate forward reads
cat "${in_dir}"/**/*"${name}"*_R1*.fastq.gz > "${out_dir}"/"${upname}"_R1.fastq.gz

## Concatenate reverse reads
cat "${in_dir}"/**/*"${name}"*_R2*.fastq.gz > "${out_dir}"/"${upname}"_R2.fastq.gz

echo
echo "End time:"
date
