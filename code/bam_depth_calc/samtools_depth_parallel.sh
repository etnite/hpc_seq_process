#!/bin/bash

## Run Samtools Depth - parallel
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script runs the samtools depth program in a single genomic region for
## a list of BAM files. It uses the associated Julia script samtools_depth_parse.jl
## to format the results. In a HPC setting, it is intended to be used with a
## parallel dispatch script to run multiple genomic regions simultaneously.
##
## The script takes two positional, command-line arguments:
##   1) The path to a .bed file defining regions of BAM files to extract.
##      This file should be sorted by chromosome and then start position, and should
##      not contain any overlapping regions.
##   2) An integer, defining which line of the .bed file to use. Supplying 0 will
##      use all lines of the .bed file. This integer will typically be supplied
##      by the parallel dispatch script.
##
## All other constants are set within the script, below this description.
##
## If using all regions, the output file will be placed in the user-defined
## output directory as "all_regions.txt". If a single region is filtered, then
## the output file will have a name formatted as [region_integer]_[chrom]_[start]_[end].txt
## The integer in the file name is formatted with leading zeroes to allow concatenation
## of multiple files without a sorting step. For instance, an output file name
## of region number 15 out of 1000 might be formatted as:
##
##   0015_2A_100000_200000.txt
##
## The script will also always output a separate header file, with region
## integer 0.
##
## The output is a pseudo .bed file with the following columns:
## 1. Chromosome
## 2. Position - 1
## 3. Position
## 4. Sum of read depths
## 5. Min read depth
## 6. Max read depth
## 7. Median read depth
## 8. Mean read depth
## 9. Read depth standard deviation
##
## The output file can then be filtered and regions with reads can be merged
## together using bedtools.
################################################################################


#### User-defined constants ####################################################

## Set parse_script to the absolute path of the companion Julia parsing script
## max_dp will cause samtools to stop counting after reaching the specified number
##   of reads in order to save time
## min_qual specifies minimum read mapping quality to include a read
## min_mdn specifies the minimum required median read depth to output a position
parse_script="/home/brian.ward/repos/hpc_seq_process/code/bam_depth_calc/samtools_depth_for_arrayer/samtools_depth_parse.jl"
bam_list_file="/home/brian.ward/samp_and_file_lists/groupA_bam_list.txt"
max_dp=200
min_qual=0
min_mdn=4
out_dir="/project/guedira_seq_map/Allegro_test/bam_depth_calc_mq0_min_mdn4"


#### Executable ################################################################


module load singularity/3.7.1
module load samtools
module load julia

echo
echo "Start time:"
date

regions_bed=$1
array_ind=$2

## This ensures output files will have names in order, to avoid needing to
## sort after concattenating them together
n=$(($(wc -l < "$regions_bed" | wc -c) - 1))
head_prefix=$(printf "%0${n}d" 0)
prefix=$(printf "%0${n}d" $array_ind)

## Create output and temp directories, only for first array value, or if array_ind is 0
if [[ $array_ind -eq 0 || $array_ind -eq 1 ]]; then
    mkdir -p "$out_dir"
    mkdir "${out_dir}/temp_files"
    echo -e "chrom\tpos0\tpos\tsum\tmin\tmax\tmedian\tmean\tsd" > "${out_dir}/${head_prefix}_header.txt"
else
    sleep 10s
fi

## If the supplied array_ind is 0, we copy full regions .bed file to temp dir.
## Otherwise extract single region from .bed file
if [[ array_ind -eq 0 ]]; then
    label="all_regions"
    cp "$regions_bed" "${out_dir}/temp_files/${label}.bed"
else
    region=$(head -n $array_ind "$regions_bed" | tail -n 1)
    suffix=$(echo "$region" | tr "\t" "_")
    label="${prefix}_${suffix}"
    echo "$region" > "${out_dir}/temp_files/${label}.bed"
fi


## Run the depth calculation and parsing
samtools depth -f "$bam_list_file" -b "${out_dir}/temp_files/${label}.bed" -d $max_dp -Q $min_qual |
    julia "$parse_script" |
    awk -v mm=$min_mdn 'BEGIN{OFS = "\t"} $7 >= mm {print $0}' > "${out_dir}/${label}.txt"

echo
echo "End time:"
date
