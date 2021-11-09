#!/bin/bash

## Merge VCF/BCF files
##
## This script will merge together a set of VCF or BCF files (or a mixture of
## the two). As always, working with BCF files will be significantly faster. It
## requires that bcftools be installed in the user's PATH.
##
## The script takes four positional input parameters:
##   1. Path to file listing input VCFs/BCFs
##   2. Path to output merged file
##   3. true/false indicating whether to only retain biallelic SNPs
##   4. Number of cores to use (for compression only)
##
## Normalization with norm_vcfs_parallel.sh may be required prior to using this
## script.
##
## NOTES ON SAMPLE NAMES
##   bcftools merge assumes that different files will contain different samples,
##   and will by default exit if the same sample is found in multiple files. In
##   this script, samples that appear more than once will have the index of their
##   file appended to their name. This can cause issues when the resulting VCF
##   file is matched up to another file, such as one containing phenotypic data.
##   You can use merge_vcf_names_resolve.R to try and fix things.
################################################################################


module load singularity/3.7.1
module load bcftools

echo
echo "Start time:"
date

## Get command line args
in_list="$1"
out_file="$2"
biallelic_snps_only="$3"
nthreads="$4"

## Get first letter of true/false options
biallelic_snps_only="${biallelic_snps_only:0:1}"

## Create indices if not present
for i in "$in_list"; do
    if [[ ! -f "${i}.csi" ]]; then bcftools index --threads $nthreads -c "$i"; fi
done

## Perform the merge
if [[ "$biallelic_snps_only" == [Tt] ]]; then
    bcftools merge --file-list "$in_list" \
        --force-samples \
        --merge both \
        --output-type u |
        bcftools view - -m2 -M2 -v snps \
            --threads $nthreads \
            --output-type b \
            --output "$out_file"
else
    bcftools merge --file-list "$in_list" \
        --threads $nthreads \
        --force-samples \
        --merge both \
        --output-type u \
        --output "$out_file"
fi
bcftools index --threads $nthreads -c "$out_file"


echo
echo "End time:"
date
