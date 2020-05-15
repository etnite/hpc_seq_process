#!/bin/bash

## Demultiplex paired fastq files by read group and convert to uBAMs
##
## This script is designed to take a pair of fastq files, and demultiplex them
## by their read groups, finally converting them to interleaved unaligned BAM
## (uBAM) files.
##
## uBAM files are  the preferred storage medium for unaligned read data in the
## GATK best practices pipeline. Specifically, there should be one uBAM per
## read group. What is a read group? Essentially, it is one run of a sequencing
## instrument. For instance, if a single sample is sequenced across multiple
## Illumina lanes, then these are separate read groups, despite being prepared
## from the same sample and library.
##
## This script is specifically for use when we have fastq files containing more
## than one read group, but we don't know what these read groups are ahead of
## time. This can happen with data downloaded from the NCBI Sequence Read Archive,
## or the European Nucleotide Archive. It can also happen when someone gives
## us single-sample fastq files that consist of multiple read-group fastq files
## concatenated together.
##
## Unfortunately in this case we must use a brute-force method to find the read
## groups, by reading every 4th line (i.e. the read ID lines) of one of the paired 
## files, parsing flowcell, lane, and barcode info, and then retaining only unique
## values. The script WILL ONLY WORK WITH 4-LINE fastqs. If you have fastqs that
## split read strings across multiple lines, that will lead to big problems
##
## BBDuk's demuxbyname.sh is used to split apart the fastq files and convert
## to interleaved BAM. Finally, picard AddOrReplaceReadGroups is used to add
## relevant read group info to the output uBAM
################################################################################



#### User-defined variables ####

fastq_dir="/project/genolabswheatphg/raw_data/v1_hapmap"
out_dir="/project/genolabswheatphg/merged_fastqs/v1_hapmap"
samp_file="/home/brian.ward/repos/wheat_phg/sample_lists/v1_hapmap_bioproj/sample_names.txt"


#### Executable ####

echo
echo "demult_fastq_by_readgroup.sh"
echo "Start time:"
date

array_ind=$1
mkdir -p "${out_dir}"

## Get sample name
samp=$(head -n "${array_ind}" "${samp_file}" | tail -n 1)

## Set forward and reverse read fastq files
fq1=$(echo "${fastq_dir}"/"${samp}"*R1.fastq.gz)
fq2=$(echo "${fastq_dir}"/"${samp}"*R2.fastq.gz)

## Convert sample name to uppercase and underscores to dashes
## Personal pref - I like samples to have dashes to make pattern recognition easier
upsamp="${samp^^}"
upsamp=$(echo "${upname}" | sed 's/_/-/g')

## Read every 4th line of the first read fastq file
## Need if statement to check for SRR files
zcat "$fq1" | sed -n '1~4p' | cut -d ":" -f 3,4,10 | uniq > fcell_lane_bar.txt

## Create comma-delim list of unique read groups
names_string=$(cat "$out_dir"/fcell_lane_bar.txt | tr '\n' ',')

## Demultiplex the fastq based on the unique flowcell, lane, barcodes
demuxbyname.sh in1="$fq1" in2="$fq2" out="$upsamp"_%.bam names="$names_string"