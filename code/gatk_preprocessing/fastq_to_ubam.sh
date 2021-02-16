#!/bin/bash

## Convert fastq files to uBAM
##
## This script is designed to take input fastq files and convert them into
## unaligned BAM (uBAM) files. 
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
## one per line. If a search pattern ends in ".fastq.gz" or "fq.gz", then we assume
## that the fastq file is either single-end or interleaved.
##
## Otherwise the script will look for "R1" and "R2" in filenames to identify
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
## output uBAM file. Metadata is collected from the first line of the first R1 fastq file,
## or the first line of the fastq file if using single-end/interleaved files.
################################################################################


#### User-defined variables ####

## Input fastq directory, output directory for .bam files, and file
## of sample names, one listed per line
fastq_dir="/bioinformatics-ward/MRAseq_v2_test_set/fastqs"
out_dir="/bioinformatics-ward/MRAseq_v2_test_set/ubams"
fq_pattern_file="/home/ward.1660/pattern_files/test_fqs_patterns.txt"


#### Executable ####

#source activate gatk4

echo
echo "fastq_to_ubam.sh"
echo "Start time:"
date

array_ind=$1
mkdir -p "${out_dir}"

## Get the pattern to search for fastq files
fq_pattern=$(head -n "${array_ind}" "${fq_pattern_file}" | tail -n 1)
#fq_pattern="TRIBUTE_sub10K_L001"
samp=$(echo "$fq_pattern" | sed 's/_.*//')
upsamp="${samp^^}"

## If search pattern ends in fastq.gz or fq.gz, then we are assuming single-ended or interleaved format
## Otherwise assume paired fastq files
if [[ "$fq_pattern" == *fastq.gz ]] || [[ "$fq_pattern" == *fq.gz ]]; then
    fq="${fastq_dir}/${fq_pattern}"
else
    fq=$(echo "${fastq_dir}/${fq_pattern}"*R1*fastq.gz "${fastq_dir}/${fq_pattern}"*R1*fq.gz)
    fq2=$(echo "${fastq_dir}/${fq_pattern}"*R2*fastq.gz "${fastq_dir}/${fq_pattern}"*R2*fq.gz)
fi

## Get flowcell, lane, and barcode metadata from R1 fastq file
## Get 10,001th line. Sometimes barcodes can contain Ns at beginning of file 
id_line=$(zcat "$fq1" | head -n 10001 | tail -n 1)
fcell=$(echo $id_line | cut -d ":" -f 3)
lane=$(echo $id_line | cut -d ":" -f 4)
bcode=$(echo $id_line | cut -d ":" -f 10)

## Create name of output .bam file
out_bam="${fq_pattern}".bam

## Run Picard FastqToSam
## Creates single-ended or interleaved BAM sorted by read names
if [[ "$fq_pattern" == *fastq.gz ]] || [[ "$fq_pattern" == *fq.gz ]]; then
    gatk FastqToSam \
        -F1 "$fq1" \
        -O "$out_dir"/"$out_bam" \
        -RG "$fcell"_"$lane" \
        -SM "$upsamp" \
        -PL ILLUMINA \
        -PU "$fcell"_"$lane"."$bcode" \
        -LB "$upsamp"_lib
else
    gatk FastqToSam \
        -F1 "$fq1" \
        -F2 "$fq2" \
        -O "$out_dir"/"$out_bam" \
        -RG "$fcell"_"$lane" \
        -SM "$upsamp" \
        -PL ILLUMINA \
        -PU "$fcell"_"$lane"."$bcode" \
        -LB "$upsamp"_lib
fi

#source deactivate
echo
echo "End time:"
date
