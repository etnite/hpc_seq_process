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
F_out_dir="/project/guedira_seq_map/Allegro_test/strand_split_bams/f_strand"
R_out_dir="/project/guedira_seq_map/Allegro_test/strand_split_bams/r_strand"
mq_thresh=0


#### Executable ################################################################

module load singularity/3.7.1
module load samtools

echo
echo "Start time:"
date

bam_list=$1
array_ind=$2

mkdir -p "$F_out_dir"
mkdir -p "$R_out_dir"

## Construct forward and reverse strand filenames
bam_file=$(head -n $array_ind "$bam_list" | tail -n 1)
bam_base=$(basename "$bam_file")
bam_base="${bam_base%.*}"
f_file="${F_out_dir}/${bam_base}_F.bam"
r_file="${R_out_dir}/${bam_base}_R.bam"

## Print info on input and output files to stdout
echo
echo "input file: $bam_file"
if [[ -f "$bam_file" ]]; then
   echo "input file exists"
else
   echo "input file does not exist"
fi

echo
echo "F out file: $f_file"
echo "R out file: $r_file"

## Split the BAM file
samtools view -h "$bam_file" -q $mq_thresh -F 20 -b -o "$f_file"
samtools view -h "$bam_file" -q $mq_thresh -f 16 -b -o "$r_file"

echo
echo "End time:"
date
