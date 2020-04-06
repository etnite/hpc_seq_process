#!/bin/bash

## Subset/Filter BAM file
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script takes a single BAM file as input, filters it, and optionally 
## subsets it based on a
## regions file. The regions file is a text file listing regions of interest,
## one per line. Valid region formats are:
##
##   1A             (All of chrom 1A)
##   1A:100000      (chrom 1A starting at 100,000bp, going to end)
##   1A:100-500     (chrom 1A from 100 to 500bp)
##
## The output consists of a single filtered, subsetted, sorted BAM file.
##
## This script is intended to be used with parallelize.sh, to enable independent
## parallel runs on multiple samples simultaneously.
##
## NOTES: 
##
##   1) This script is intended to be modified by hand as needed, so some
##      flags/arguments in the samtools view and samtools index calls should be
##      inspected and changed as necessary
##   2) Input BAM files should be sorted
##   3) Working directory inherited from parallelizing script - it is easiest
##      to define absolute paths
################################################################################


#### User-defined constants ####

in_dir=""
out_dir=""
regions_file=""


#### Executable  ####

module load samtools

echo
echo "Start subset_filt_bam.sh"
echo "Start time:"
date

mkdir -p "${out_dir}"
samp=$1

## Input and output BAM paths
in_bam="${in_dir}"/"${samp}".bam
out_bam="${out_dir}"/"${samp}".bam

## Customize as necessary
samtools view "${in_bam}" -h -R "${regions_file}" -f 2 -b -o "${out_bam}"
samtools index -c "${out_bam}"

echo
echo "End time:"
date
