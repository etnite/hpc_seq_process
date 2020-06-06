#!/bin/bash
set -e


## Parallel mpileup and variant calling
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script performs variant calling using BCFTools, running the bcftools
## mpileup and call commands for a single region. Note that this script DOES NOT
## use .bed files to define regions. Rather, it uses samtools/bcftools region
## files, where regions are defined as:
##
##    chromosome:from-to
##
## One region per line. Unlike .bed files, both the from and to positions are
## 1-based. 
##
## This script performs variant calling on a single region. It is designed to be 
## used with a parallel dispatching script, such as code/parallel_dispatch/arrayer.sh See
## code/make_regions_file.py for a script that creates non-overlapping windows
## across a reference genome. 
##
## A minimum mapping-quality threshold may be set. For bowtie2, probability of
## an incorrect alignment (p) is related to the mapping quality value (Q) by:
##
##   p = 10 ^ -(Q/10)
##
## A mapq score of 5 translates to a ~32% chance of misalignment
## A mapq score of 10 translates to a 10% chance
## A mapq score of 20 translates to a 1% chance, etc.
##
## Currently the script only uses reads mapped in "proper pairs". This can be disabled
## by deleting the "--rf 2" in the bcftools mpileup call.
################################################################################


#### User-defined constants ####

bams_list="/home/brian.ward/US_excap_bam_list.txt"
ref_gen="/project/genolabswheatphg/v1_refseq/whole_chroms/Triticum_aestivum.IWGSC.dna.toplevel.fa"
regions_file="/home/brian.ward/region_files/v1_100Mb_window_regions.txt"
out_dir="/project/guedira_seq_map/brian/US_excap/v1_variants/region_bcfs"
mq_val=20


#### Executable ####

module load bcftools

echo
echo "Start time:"
date

mkdir -p "$out_dir"
array_ind=$1


## Get region
region=$(head -n "$array_ind" "$regions_file" | tail -n 1)
echo
echo "Input genomic region: ${region}"

## Reformat region for output file
out_reg=$(echo "$region" | sed 's/:/_/' | sed 's/-/_/')

## Perform the mpileup and calling
bcftools mpileup --fasta-ref "$ref_gen" \
                 --bam-list "$bams_list" \
                 --min-MQ "$mq_val" \
                 --regions "$region" \
                 --incl-flags 2 \
                 --annotate FORMAT/DP,FORMAT/AD \
                 --output-type u |
   bcftools call --multiallelic-caller \
                 --variants-only \
                 --output-type b \
                 --output "${out_dir}/${out_reg}.bcf"

echo
echo "End time:"
date
