#!/bin/bash

## Convert paired fastq files to uBAMs
##
## This script is designed to take a pair of fastq files, and convert them into a
## single interleaved unaligned BAM (uBAM) file.
##
## uBAM files are  the preferred storage medium for unaligned read data in the
## GATK best practices pipeline. Specifically, there should be one uBAM per
## read group. What is a read group? Essentially, it is one run of a sequencing
## instrument. For instance, if a single sample is sequenced across multiple
## Illumina lanes, then these are separate read groups, despite being prepared
## from the same sample and library.
##
## If you have fastq files containing multiple read groups, then make sure to use
## demux_fastq_by_readgroup.sh instead.
##
## The script assumes that a read group has a single pair of fastq files. The user
## supplies a text file which lists the patterns in fastq file names to look for,
## one per line. The script will look for "R1" and "R2" in filenames to identify
## the read pairs. For instance, a typical entry might be:
##
## Zenda_S29_L001
##
## Which would match the paired fastq files:
##
## Zenda_S29_L001_R1_001.fastq.gz
## Zenda_S29_L001_R2_001.fastq.gz
##
## The script will take everything up to the first underscore as the sample name.
## The sample name is converted to uppercase and stored in the metadata of the
## output uBAM file. Metadata collected from the first line of the first R1 fastq file.
################################################################################



#### User-defined variables ####

## Input fastq directory, output directory for .bam files, and file
## of sample names, one listed per line
fastq_dir="/project/genolabswheatphg/merged_fastqs/SRW_excap"
out_dir="/project/genolabswheatphg/merged_fastqs/SRW_excap/fq_demux_test"
fq_pattern_file="/home/brian.ward/repos/wheat_phg/sample_lists/v1_hapmap_bioproj/sample_names.txt"


#### Executable ####

source activate gatk4

echo
echo "fastq_to_ubam.sh"
echo "Start time:"
date

array_ind=$1
mkdir -p "${out_dir}"

## Get the pattern to search for fastq files
#fq_pattern=$(head -n "${array_ind}" "${fq_pattern_file}" | tail -n 1)
fq_pattern="TRIBUTE_sub10K_L001"
samp=$(echo "$fq_pattern" | sed 's/_.*//')
upsamp="${samp^^}"

## Set forward and reverse read fastq files
fq1=$(echo "${fastq_dir}"/"${fq_pattern}"*R1*fastq.gz)
fq2=$(echo "${fastq_dir}"/"${fq_pattern}"*R2*fastq.gz)

## Get flowcell, lane, and barcode metadata from R1 fastq file
## Get 10,001th line. Sometimes barcodes can contain Ns at beginning of file 
id_line=$(zcat "$fq1" | head -n 10001 | tail -n 1)
fcell=$(echo $id_line | cut -d ":" -f 3)
lane=$(echo $id_line | cut -d ":" -f 4)
bcode=$(echo $id_line | cut -d ":" -f 10)

## Create name of output .bam file
out_bam="${fq_pattern}".bam

## Run Picard FastqToSam
## Creates interleaved BAM sorted by read names
gatk FastqToSam \
    -F1 "$fq1" \
    -F2 "$fq2" \
    -O "$out_dir"/"$out_bam" \
    -RG "$fcell"_"$lane" \
    -SM "$upsamp" \
    -PL ILLUMINA \
    -PU "$fcell"_"$lane"."$bcode" \
    -LB "$upsamp"_lib

source deactivate
echo
echo "End time:"
date
