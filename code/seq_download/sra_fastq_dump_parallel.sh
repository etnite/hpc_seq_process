#!/bin/bash

## Dump SRA files to FASTQ
##
## This script takes two inputs:
##   1) Two-column tab-delimited file, with SRR numbers in first column (e.g.
##      "SRR1211028", and sample name in second column. This file would typically
##      be produced using get_sra_nums.sh
##   2) Path to directory to write output
##
## This script is designed to be used with arrayer.sh, which passes a single
## integer as a positional argument. This integer is used to select a single
## line from the user-defined two-column sample file.
################################################################################


#### User-defined constants ####

samp_file="/home/brian.ward/repos/wheat_phg/sample_lists/v1_hapmap_bioproj/SRA_samp_list.tsv"
out_dir="/project/genolabswheatphg/raw_data/v1_hapmap"


#### Executable ####

module load sratoolkit

echo
echo "Start sra_fastq_dump_parallel.sh"
echo "Start time:"
date

mkdir -p "${out_dir}"

## Get the SRR number
srr=$(head -n "${1}" "${samp_file}" | tail -n 1 | cut -f 1)

## Get sample name
samp=$(head -n "${1}" "${samp_file}" | tail -n 1 | cut -f 2)

## Echo SRR, sample to stdout and stderr
echo
echo >&2
echo "SRR number: ${srr}"
echo "SRR number: ${srr}" >&2
echo "Sample name: ${samp}"
echo "Sample name: ${samp}" >&2

## Get fastq files
fastq-dump --gzip \
           --skip-technical \
           --readids \
           --read-filter pass \
           --dumpbase \
           --split-3 \
           --outdir "${out_dir}" \
           --clip "${srr}"

## Rename fastq files to include sample name
mv "${out_dir}"/"${srr}"_pass_1.fastq.gz "${out_dir}"/"${samp}"_R1.fastq.gz
mv "${out_dir}"/"${srr}"_pass_2.fastq.gz "${out_dir}"/"${samp}"_R2.fastq.gz

echo
echo "End time:"
date

