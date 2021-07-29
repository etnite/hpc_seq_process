#!/bin/bash

## Split BAM file by strand
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script filters a BAM (or SAM) file, splitting alignments to the positive
## and negative strands into separate BAM files. The two strand files are
## placed in separate output directories. Note that unaligned reads are not
## included in either output file.
##
## The script takes two positional, command-line arguments:
##   1) The path to a .txt file listing input BAM files, one per line
##   2) An integer, defining which line of the BAM list to use. Supplying 0 will
##      use all lines of the file. This integer will typically be supplied
##      by a parallel dispatch script, which will allow for filtering multiple
##      regions simultaneously.
################################################################################


#### User-Supplied Constants ###################################################

## Forward and reverse strand output directories
## Optional minimum read quality threshold - set to 0 to disable
F_out_dir=""
R_out_dir=""
mq_thresh=20


#### Executable ################################################################

echo
echo "Start time:"
date

bam_list=$1
array_ind=$2

## Construct forward and reverse strand filenames
bam_file=$(head -n $arr_ind "$bam_list" | tail -n 1)
bam_base=$(basename "$bam_file")
bam_base="${bam_base%.*}"
f_file="${F_out_dir}/${bam_base}_F.bam"
r_file="${R_out_dir}/${bam_base}_R.bam"

## Split the BAM file
samtools view -h "$bam_file" -q $mq_thresh -F 20 -b -o "$f_file"
samtools view -h "$bam_file" -q $mq_thresh -f 16 -b -o "$r_file"

echo
echo "End time:"
date
