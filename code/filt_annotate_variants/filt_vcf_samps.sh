#!/bin/bash

## Filter a VCF/BCF file sample-wise
##
## This script will filter a VCF or BCF file in a sample-wise fashion. Specifically,
## it implements the following filters, all optional:
##   1. Retain only samples in a list supplied by the user
##   2. Retain samples with missingness below a user-defined threshold
##   3. Retain samples with heterozygosity below a user-defined threshold
##
## It takes the path to the input file as a positional argument. The output file
## is placed in the same directory with the input file, with "samp_filt" inserted
## before the extension. The file will output either a bgzipped VCF file or a BCF
## file depending on the type of file supplied as input.
##
## NOTE that the filtering will be slightly inaccurate in the presence of indels,
## as the bcftools stats command reports statistics on these separately from
## the SNP per-sample counts that are utilized in this script.
###############################################################################


#### User-defined Constants ############

samp_file="none"
max_miss=0.85
max_het=0.3


#### Executable ########################

echo
echo "Start time:"
date

## Path handling
vcf_in=$1
out_dir=$(dirname "$vcf_in")
base="${vcf_in%%.*}"
base="${vcf_in%%.*}"

ext="${vcf_in#*.}"
if [[ "$ext" == "gz" ]]; then ext="vcf.gz"; fi
vcf_out="${base}_samp_filt.${ext}"

## Write filtering parameters to output directory
echo -e "Input VCF\t${vcf_in}" > "${base}_samplewise_filt_params.txt"
echo -e "Output VCF\t${vcf_out}" >> "${base}_samplewise_filt_params.txt"
echo -e "Sample subset list\t${samp_file}" >> "${base}_samplewise_filt_params.txt"
echo -e "Max. missing threshold\t${max_miss}" >> "${base}_samplewise_filt_params.txt"
echo -e "Max. het. threshold\t${max_het}" >> "${base}_samplewise_filt_params.txt"

## Set the file type for output file
if [[ "$ext" == "bcf" ]]; then
    ftype="b"
elif [[ "$ext" == "vcf.gz" || "$ext" == "vcf" ]]; then
    ftype="z"
fi

## Create a temp directory within the working directory
temp_dir="$(mktemp -d -p "$out_dir")"
if [[ ! "$temp_dir" || ! -d "$temp_dir" ]]; then
    echo "Could not create temporary directory"
    exit 1;
fi

## Add leading/trailing tabs to sample names
if [[ -f "$samp_file" ]]; then
    awk 'BEGIN{OFS = ""} {print "\t", $1, "\t"}' "$samp_file" | sort > "${temp_dir}/samp_file.txt"
else
    bcftools query --list-samples "$vcf_in" | awk 'BEGIN{OFS = ""} {print "\t", $1, "\t"}' | sort > "${temp_dir}/samp_file.txt"
fi

## Calculate statistics from input VCF/BCF
bcftools stats -s - "$vcf_in" > "${temp_dir}/stats.txt"

## Filtering block - isolate per-sample counts (PSC) from stats file
## Then retain samples in user-defined file
## Then retain samples passing max_miss threshold
## Finally retain samples passing max_het threhsold
echo
echo "Filtering VCF..."
grep "^PSC" "${temp_dir}/stats.txt" |
    fgrep -f "${temp_dir}/samp_file.txt" |
    awk -v miss=$max_miss '$14/($4 + $5 + $6 + $14) <= miss {print $0}' |
    awk -v het=$max_het '$6/($4 + $5 + $6) <= het {print $3}' > "${temp_dir}/retain.txt"

## Generate output file, only retaining passing samples
bcftools view "$vcf_in" -S "${temp_dir}/retain.txt" -O "$ftype" -o "$vcf_out"
bcftools index "$vcf_out"

## Remove temp. directory
rm -rf "$temp_dir"
echo
echo "Removed temp directory $temp_dir"

echo
echo "End time:"
date
