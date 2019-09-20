#!/bin/bash
set -e

## Filter BAM files created by Bowtie2
##
## This script is designed to work with arrayer.sh - it takes a single command
## line argument, which is an integer. It uses this integer to select a single
## .bam file from in_dir, then filters the .bam file and writes the output to
## out_dir. 
##
## Notes on isolating "uniquely mapping" reads
##
## This is not a straightforward topic, especially because the concept of
## "uniquely mapping" is ill-defined. Bowtie2 reports a MAPQ score (Q) for each
## alignment which is defined by the probability (p) that an alignment does not
## correspond to the read's true point of origin:
##
##   Q = -10 * log10(p)
##   or:
##   p = 10 ^ -(Q/10)
##
## The value of p is calculated based on the number of alternative alignments,
## and how well these other alignments match the reference. Note that in most
## cases bowtie2 DOES NOT perform an exhaustive search of alignments for each
## read. The user can specify an exhaustive search, but this is ill-advised for
## large genomes (I'm guessing anything larger than bacterial genomes). Therefore,
## the mapq value is based on a sample of the larger set of all possible
## alignments for a given read.
##
## A mapq score of 5 translates to a ~32% chance of misalignment
## A mapq score of 10 translates to a 10% chance
## A mapq score of 20 translates to a 1% chance, etc.
################################################################################ 


#### User-Supplied Constants ####

#in_dir="/project/genolabswheatphg/alignments/SRW_wholechrom_bw2_bams"
#out_dir="/project/genolabswheatphg/alignments/SRW_filt_bams"
#mq_thresh=20

in_dir="/project/genolabswheatphg/alignments/KS_HRW_wholechrom_bw2_bams"
out_dir="/project/genolabswheatphg/alignments/KS_HRW_filt_bams"
mq_thresh=20


#### Executable ####

module load samtools

mkdir -p "${out_dir}"

## Read in array index (integer) - get corresponding .bam file name
arr_ind=$1
bam_file=$(ls -1 "${in_dir}"/*.bam | head -n $arr_ind | tail -n 1)
bam_base=$(basename "${bam_file}")

## Filter input BAM file, sort and output
samtools view -h "${bam_file}" \
    -q $mq_thresh \
    -f 2 |
    samtools sort -O BAM - -o "${out_dir}"/"${bam_base}"

## Index output BAM
samtools index -c "${out_dir}"/"${bam_base}"
