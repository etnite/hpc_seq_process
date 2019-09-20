#!/bin/bash

## Subset/Downsample FASTQ file
##
## This script takes a single BAM file as input, and uses the read names in this
## file to subset a pair of mated fastq files (which were probably used to create
## the BAM file in the first place. This is primarily for testing/simulation
## purposes. It requires the both bbtools, and the program seqtk, installed in
## a conda environment called "seqtk"
##
## The script takes a string as a single positional input (NAME). It then expects to
## find:
##
##   Template bam file located at bam_dir/NAME.bam
##   Corresponding FASTQ files located at fastq_dir/NAME_R[12].fastq.gz
##
## The output consists of a single pair of subsetted fastq files, and possible
## a pair of downsampled fastq files.
##
## This script is intended to be used with parallelize.sh, to enable independent
## parallel runs on multiple samples simultaneously.
##
## NOTES:
##   1) If downsamp_dir is set to "NA", no subsetting will be performed
##   2) Working directory inherited from parallelizing script - it is easiest
##      to define absolute paths
################################################################################


#### User-defined constants ####

fastq_dir=""
bam_dir=""
subset_dir=""
downsamp_dir="NA"
downsamp_n=1000000


#### Executable  ####

module load samtools
module load bbtools
if [[ "${downsamp_dir}" != "NA" ]]; then source activate seqtk; fi

echo
echo "Start subset_downsamp_fastq.sh"
echo "Start time:"
date

samp=$1

## Isolate only properly paired reads
## Not sure if bbtools can use a BAM file or not
samtools view -h -f 2 "${bam_dir}"/"${samp}".bam > "${bam_dir}"/"${samp}"_temp.sam

## Perform the FASTQ subsetting
filterbyname.sh in="${fastq_dir}"/"${samp}"_R1.fastq.gz \
                in2="${fastq_dir}"/"${samp}"_R1.fastq.gz \
                out="${subset_dir}"/"${samp}"_R1.fastq.gz \
                out2="${subset_dir}"/"${samp}"_R2.fastq.gz \
                names="${bam_dir}"/"${samp}"_temp.sam \
                include=t

## Optionally downsample the subset fastq files
if [[ "${downsamp_dir}" != "NA" ]]; then
    seqtk sample -s10 "${subset_dir}"/"${samp}"_R1.fastq.gz $downsamp_n > "${subset_dir}"/"${samp}"_${downsamp_n}_R1.fastq.gz
    seqtk sample -s10 "${subset_dir}"/"${samp}"_R2.fastq.gz $downsamp_n > "${subset_dir}"/"${samp}"_${downsamp_n}_R2.fastq.gz
fi

## Cleanup and exit
rm "${bam_dir}"/"${samp}"_temp.sam
source deactivate

echo
echo "End time:"
date
