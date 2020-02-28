#!/bin/bash

## Convert BAM files back to fastqs
##
## The purpose of this script is to take a BAM file, and convert it back into
## paired FASTQ files. Why would you want to do this? Mainly for testing
## purposes to compare different alignment techniques.
##
## It requires that the user supply a text file with the .bam files of interest
## listed, one per line. This file would typically contain absolute paths, though
## that need not necessarily be the case. Corresponding paired fastq files are
## then output into a specified folder
##
## The samtools view call can be hand-tweaked to implement any desired filtering.
################################################################################      


bams_file="/home/brian.ward/repos/hpc_seq_process/sample_lists/SRW_bam_list.txt"
out_fastqs_dir="/project/genolabswheatphg/test_1A_fastqs"


#### Executable ####

echo
echo "Start bam_to_fastq_parallel.sh"
echo "Start time:"
date

module use /shared/7/modulefiles
module load samtools
samtools --version

array_ind=$1
mkdir -p "${out_fastqs_dir}"

## Get sample name from .bam file
## Assumes file is named "<sample_name>.bam"
bam=$(head -n "${array_ind}" "${bams_file}" | tail -n 1)
samp=$(basename "$bam")
samp="${samp%.*}"

## Set names of output fastq files
fq1="$out_fastqs_dir"/"$samp"_R1.fq.gz
fq2="$out_fastqs_dir"/"$samp"_R2.fq.gz

## Convert .bam file to .fastqs
## This call to samtools view should be customized to subset reads as desired
samtools view -q 20 "$in_bams_dir"/"$bam" 1A |
	samtools fastq - -1 "$fq1" -2 "$fq2" -0 /dev/null -s /dev/null

echo
echo "End time:"
date
