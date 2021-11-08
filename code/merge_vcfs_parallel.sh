#!/bin/bash

## Merge together multiple VCF/BCF files
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script will merge together two or more VCF or BCF files (or a mix of
## VCF and BCF files). This is a "horizontal" merge which is intended to handle
## files containing different sets of samples. For a "vertical" merge (i.e.
## files containing different sets of variants for the same samples) something
## like bcftools concat is typically used.
##
## NOTES ON SAMPLE NAMES
##   bcftools merge assumes that different files will contain different samples,
##   and will by default exit if the same sample is found in multiple files. In
##   this script, samples that appear more than once will have the index of their
##   file appended to their name. This can cause issues when the resulting VCF
##   file is matched up to another file, such as one containing phenotypic data.
##   You can use merge_vcf_names_resolve.R to try and fix things.
##

################################################################################


##### User-Defined Constants ##########

out_dir=""
normalize="true"
ref_genome=""
biallelic_snps_only="true"
keep_normed="true"


##### Executable #######################

module load singularity/3.7.1
module load bcftools

echo
echo "Start time:"
date

in_list="$1"
array_ind="$2"

## Get first letter of true/false normalize option and
normalize="${normalize:0:1}"

## Set up the output directories
mkdir -p "$out_dir"
if [[ "$normalize" == [Tt] ]]; then
	mkdir "${out_dir}/normalized_vcfs"
fi
if [[ $array_ind -eq 0 || $array_ind -eq 1 ]]; then
    rm -f "${out_dir}/merge_list.txt"
else
    sleep 20s
fi

## Construct an array - either the entire input list, or else one line of it
if [[ $array_ind -eq 0 ]]; then
    mapfile -t vcf_array < "$in_list"
else
    n_toss=(( $array_ind - 1 ))
    mapfile -t -n $array_ind -s $n_toss vcf_array < "$in_list"   ## Should work for -s 0
fi

## Now loop through our array of input files (possibly just one)
for i in "${vcf_array[@]}"; do
    if [[ ! -f "${i}.csi" ]]; then bcftools index -c "$i"; fi

    if [[ "$normalize" == [Tt] ]]; then

        ## Get the basename of the input VCF/BCF
        ## Then we create a corresponding filename for after normalization
        ## This will always be a .bcf file and may (possibly) need its extension changed
        base=$(basename "$i")
        normed_file="${out_dir}/normalized_vcfs/${base}"
        normed_file="${normed_file%.*}"
        normed_file="${normed_file%.*}"
        normed_file="${normed_file}.bcf"
        echo "$normed_file" >> "${out_dir}/merge_list.txt"

        bcftools norm "$i" \
            --check-ref s \
            --multiallelics +both \
            --output-type b \
            --output "$normed_file"
        bcftools index -c "$normed_file"

    else
        echo "$i" >> "${out_dir}/merge_list.txt"
    fi
done

sort "${out_dir}/merge_list.txt" -o "${out_dir}/merge_list.txt"

## Perform the merge
if [[ "$biallelic_snps_only" == [Tt] ]]; then
    bcftools merge --file-list "${out_dir}/merge_list.txt" \
        --force-samples \
        --merge both \
        --output-type u |
        bcftools view - -m2 -M2 -v snps \
            --output-type u \
            --output "${out_dir}/merged.bcf"
else
    bcftools merge --file-list "${out_dir}/merge_list.txt" \
        --force-samples \
        --merge both \
        --output-type u \
        --output "${out_dir}/merged.bcf"
fi
bcftools index -c "${out_dir}/merged.bcf"

if [[ ! "$keep_normed" == [Tt] ]]; then rm -rf "${out_dir}/normalized_vcfs"; fi

echo
echo "End time:"
date
