#!/bin/bash

## Download FASTQ files from ENA
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script is designed to download paired FASTQ files from the ENA using an
## FTP connection (though it can also use Aspera if configured). This should be
## faster than downloading from the NCBI sequence read archive, as it avoids
## having to convert SRA to FASTQ format. In addition, FASTQ files from the ENA
## should be properly formatted and paired. In addition to FASTQ files, this 
## script will download an .xml file containing meta info for each sample
##
## The script takes two inputs:
##   1) Two-column tab-delimited file, with SRR numbers in first column (e.g.
##      "SRR1211028", and sample name in second column. This file would typically
##      be produced using get_sra_nums.sh
##   2) Path to directory to write output - note that this script will produce
##      a subdirectory for each sample within the output directory
##
## This script is designed to be used with arrayer.sh, which passes a single
## integer as a positional argument. This integer is used to select a single
## line from the user-defined two-column sample file.
##
## This script relies on a Bioconda environment called "ena_dl" which contains
## the package "enabrowsertools". To create this, run:
##
##   conda create --name ena_dl enabrowsertools
## 
## and follow prompts
################################################################################

#### User-defined constants ####

samp_file="/home/brian.ward/repos/wheat_phg/sample_lists/v1_hapmap_bioproj/SRA_repeats.tsv"
out_dir="/project/genolabswheatphg/raw_data/v1_hapmap"


#### Executable ####

module load miniconda
source activate ena_dl  ## Newer conda versions use conda activate

echo
echo "Start ena_fastq_dl.sh"
echo "Start time:"
date

mkdir -p "${out_dir}"

## Get the SRR number
srr=$(head -n "${1}" "${samp_file}" | tail -n 1 | cut -f 1)

## Get sample name
samp=$(head -n "${1}" "${samp_file}" | tail -n 1 | cut -f 2)

## Download fastqs and metadata
enaDataGet --meta --format fastq --dest "${out_dir}" "${srr}"

## Rename fastq files to include sample name
mv "${out_dir}"/"${srr}"/"${srr}"_1.fastq.gz "${out_dir}"/"${srr}"/"${samp}"_R1.fastq.gz
mv "${out_dir}"/"${srr}"/"${srr}"_2.fastq.gz "${out_dir}"/"${srr}"/"${samp}"_R2.fastq.gz

source deactivate  ## Newer conda versions use conda deactivate

echo
echo "End time:"
date
