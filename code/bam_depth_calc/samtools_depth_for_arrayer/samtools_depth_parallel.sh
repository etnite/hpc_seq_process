#!/bin/bash

## Run Samtools Depth - parallel
################################################################################


#### User-defined constants ####################################################

## Set parse_script to the absolute path of the companion Julia parsing script
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

rm -rf "${out_dir}/temp_files"

echo
echo "End time:"
date
