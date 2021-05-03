#!/bin/bash

## Merge BAM files from same sample
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script will recursively search for alignment (BAM) files within an input
## directory, and will then merge together all files whose names contain a match 
## to a pattern string read in from a file (one pattern per line).
##
## This script is meant to be used in conjunction with arrayer.sh to enable
## parallel execution on multiple samples simultaneously. arrayer.sh will then
## handle the parallel dispatch of this script for each line in the input file
## specifying sample names.
##
## NOTES: 
##   1) Generally, the search patterns should be sample names. The script assumes
##      that the sample names occur at the beginning of the basename of each file,
##      up to the first underscore. e.g. for the file:
##         
##          input_directory/subdirectory/name_S4_L001.bam
##      
##      "name" is the sample name. The script prepends "/" and appends "_" to the
##      search pattern to try and increase specificity.
##   2) The user should be sure that they actually want these alignments merged,
##      since the script will do so without any regard for the contents of the
##      files. For example, it will merge alignments from different read groups
##      (i.e. flowcell/lane combo) and alignments made with different aligners.
##   3) The script will overwrite any preexisting .bam file located at the output
##      path.
################################################################################      


#### User-Defined Constants ####

in_bams_dir="/project/genolabswheatphg/US_excap/v1_alignments"
out_bams_dir="/project/genolabswheatphg/US_excap/v1_merged_mq20_bams"
patterns_file="/home/brian.ward/search_pattern_files/ags_la_samp_names.txt"

## Convert sample name in RG line of BAM header to uppercase? (TRUE/FALSE)
name2upper="TRUE"

## Minimum mapping quality to use for filtering reads
## Set to 0 to disable any filtering
min_qual=20


#### Executable ####

echo
echo "Start merge_bams_parallel.sh"
echo "Start time:"
date

module load samtools

## Input sanity check
name2upper="${name2upper^^}"
if [[ "$name2upper" != "TRUE"  ]] && [[ "$name2upper" != "FALSE" ]]; then
    echo
    echo "Error - please set name2upper to either 'TRUE' or 'FALSE'"
    echo "It is currently set to: ${name2upper}"
    exit 1;
fi

array_ind=$1
mkdir -p "$out_bams_dir"

## Get the input pattern
patt=$(head -n "$array_ind" "$patterns_file" | tail -n 1)
echo
echo "Input search pattern: ${patt}"

## Prepend "/" and "_" to pattern if necessary for greater specificity
if [[ "$patt" != /* ]]; then patt="/${patt}"; fi
if [[ "$patt" != *_ ]]; then patt="${patt}_"; fi

## Strip out the "/" and "_" to isolate sample name
samp=$(echo "$patt" | sed 's/^\///' | sed 's/_$//')

## Convert sample name to uppercase if specified
if [[ "$name2upper" == "TRUE" ]]; then samp="${samp^^}"; fi

## Recursively find all .bam files in the input directory
in_bams=( $(find "$in_bams_dir" -name '*.bam') )

## Dump names of bams matching sample name pattern into a temporary file
printf '%s\n' "${in_bams[@]}" | grep "$patt" > "${out_bams_dir}/${samp}_bam_list.txt"

## Merge and sort the BAM files
## This step currently does not perform any filtering, but can be customized
## with a call to samtools view
samtools merge -f -b "${out_bams_dir}/${samp}_bam_list.txt" "${out_bams_dir}/${samp}.bam"
	
## Filter by mapping quality if threshold is set to greater than 0
if [[ "$min_qual" -gt 0 ]]; then
    samtools view -b -q "$min_qual" "${out_bams_dir}/${samp}.bam" > "${out_bams_dir}/${samp}_temp_minqual.bam"
    mv "${out_bams_dir}/${samp}_temp_minqual.bam" "${out_bams_dir}/${samp}.bam"
fi

## Index the merged bam file
samtools index -c "${out_bams_dir}/${samp}.bam"

rm "${out_bams_dir}/${samp}_bam_list.txt"


echo
echo "End time:"
date
