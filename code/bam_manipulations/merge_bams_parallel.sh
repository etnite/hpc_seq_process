#!/bin/bash

## Merge BAM files from same sample
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script will recursively search for alignment (BAM) files within an input
## directory, and will then merge together all files whose names contain a match 
## to a sample name read in from a file (one sample name per line).
##
## This script is meant to be used in conjunction with arrayer.sh to enable
## parallel execution on multiple samples simultaneously. arrayer.sh will then
## handle the parallel dispatch of this script for each line in the input file
## specifying sample names.
##
## NOTES: 
##   1) The user should take care to avoid "collisions" between sample names.
##      For instance, specifying different samples "Virginia14" and "Virginia14-01"
##      would cause problems.
##   2) This script assumes that the various BAM files being merged have the same
##      read group defined in their headers. As long as this is really the case,
##      it will just use the common read group, without appending any unique
##      suffixes. However, if different files being merged have different read
##      groups, these will be preserved (this is not a bug, but probably not
##      what the user wants in this case).
##   3) Generalizing further, the user should be certain that these are alignments
##      that should actually be merged - i.e. from the same sample, using the
##      same aligner, etc.
################################################################################      


in_bams_dir="/project/genolabswheatphg/alignments/ERSGGL_SRW_bw2_bams"
out_bams_dir="/project/genolabswheatphg/alignments/ERSGGL_SRW_bw2_merged_excap_GBS_bams"
samples_file="/home/brian.ward/repos/wheat_phg/sample_lists/SS-MPV57.txt"


#### Executable ####

echo
echo "Start merge_bams_parallel.sh"
echo "Start time:"
date

module load samtools
samtools --version

array_ind=$1
mkdir -p "${out_bams_dir}"

## Get sample name
samp=$(head -n "${array_ind}" "${samples_file}" | tail -n 1)

## Recursively find all .bam files in the input directory
shopt -s globstar nullglob
in_bams=( "${in_bams_dir}"/**/*.bam )

## Dump names of bams matching sample name pattern into a temporary file
printf '%s\n' "${in_bams[@]}" | grep "${samp}" > "${out_bams_dir}"/"${samp}"_bam_list.txt

## Merge and sort the BAM files
## This step currently does not perform any filtering, but can be customized
## with a call to samtools view
samtools merge -c -b "${out_bams_dir}"/"${samp}"_bam_list.txt "${out_bams_dir}"/"${samp}".bam
	
## Index the merged bam file
samtools index -c "${out_bams_dir}"/"${samp}".bam

rm "${out_bams_dir}"/"${samp}"_bam_list.txt


echo
echo "End time:"
date
