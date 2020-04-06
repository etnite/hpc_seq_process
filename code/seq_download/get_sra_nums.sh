#!/bin/bash

## Get SRA samples associated with a BioProject
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script accepts two positional input parameters:
##   1) The accession number of an SRA bioproject
##   2) Output directory
##
## The outputted list may need some additional formatting. In addition, a test
## with the v1 wheat hapmap (bioproject SRP032974) showed that some samples had
## multiple associated SRA files - right now these are just given a unique
## counter number to keep files from being overwritten inadvertently.
################################################################################


#module load edirect

mkdir -p "${2}"

## Get SRA run info (.csv file with lots of columns)
esearch -db sra -query "${1}" |
    efetch --format runinfo > "${2}"/SRA_runinfo.csv

## Isolate and format SRR numbers and sample names
cut -d "," -f 1,30 "${2}"/SRA_runinfo.csv |
    grep -v "Run" |
    sed 's/SeqCap_//' |
    sed '/^[[:space:]]*$/d' |
    sort -k1,1 -k2,2 |
    uniq |
    tr '[:lower:]' '[:upper:]' |
    sed 's/_/-/g' |
    tr ',' '\t' > "${2}"/SRA_samp_list.tsv

## Get list of sample names
cut -f 2 "${2}"/SRA_samp_list.tsv | uniq > "${2}"/sample_names.txt

## Append line number to sample names to make them unique
## (Avoids possible overwriting of files)
awk '{print $s "_samp" NR}' "${2}"/SRA_samp_list.tsv > "${2}"/temp_samp_list.tsv
mv "${2}"/temp_samp_list.tsv "${2}"/SRA_samp_list.tsv
