#!/bin/bash
set -e

## Run Samtools Depth - parallel
################################################################################


#### User-defined constants ####################################################

## Note that SNP depth and proportion of missing data are highly correlated

parse_script=""
bam_list_file=""
max_dp=200
min_qual=0
out_dir=""


#### Executable ################################################################



echo
echo "Start time:"
date

regions_bed=$1
array_ind=$2

## Create output and temp directories, only for first array value, or if array_ind is 0
if [[ $array_ind -eq 0 || $array_ind -eq 1 ]]; then
    mkdir -p "$out_dir"
    mkdir "${out_dir}/temp_files"
else
    sleep 10s
fi

## If the supplied array_ind is 0, we copy full regions .bed file to temp dir.
## Otherwise extract single region from .bed file
if [[ array_ind -eq 0 ]]; then
    label="all_regions"
    cp "$regions_bed" "${out_dir}/temp_files/${label}.bed"
else
    ## This ensures output files will have names in order, to avoid needing to
    ## sort after concattenating them together
    n=$(($(wc -l < "$regions_bed" | wc -c) - 1))
    prefix=$(printf "%0${n}d" $array_ind)

    region=$(head -n $array_ind "$regions_bed" | tail -n 1)
    suffix=$(echo "$region" | tr "\t" "_")
    label="${prefix}_${suffix}"
    echo "$region" > "${out_dir}/temp_files/${label}.bed"
fi


samtools depth -f "$bam_list_file" -b "${out_dir}/temp_files/${label}.bed" -d $max_dp -Q $min_qual |
    "$parse_script" > "${out_dir}/${label}.txt"


rm -rf "${out_dir}/temp_files"

echo
echo "End time:"
date
