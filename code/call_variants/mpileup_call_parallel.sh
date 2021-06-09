#!/bin/bash

## Parallel mpileup and variant calling
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script performs variant calling using BCFTools,
##
## It takes two positional arguments. The first is a .bed file of genomic regions
## and the second is an integer selecting which region to
## perform mpileup/variant calling on. If 0 is provided as the integer, then
## variant calling is performed on all regions in the file.
##
## The script is designed to be used with a parallel dispatching script, such as
## code/parallel_dispatch/arrayer.sh See code/make_regions_bed.py for a script
## that creates non-overlapping windows across a reference genome.
##
## A minimum mapping-quality threshold (mq) may be set. For bowtie2, probability of
## an incorrect alignment (p) is related to the mapping quality value (Q) by:
##
##   p = 10 ^ -(Q/10)
##
## A mapq score of 5 translates to a ~32% chance of misalignment
## A mapq score of 10 translates to a 10% chance
## A mapq score of 20 translates to a 1% chance, etc.
##
## Currently the script only uses reads mapped in "proper pairs". If using
## single-ended sequencing this will remove all reads. This behavior can be disabled
## by deleting the "--incl-flags 2" in the bcftools mpileup call.
################################################################################


#### User-defined constants ####

bams_list="/home/brian.ward/samp_and_file_lists/groupA_bam_list.txt"
ref_gen="/project/guedira_seq_map/ref_genomes/v1_refseq_w_KIMs/CSv1_refseq_w_KIMs.fa"
out_dir="/project/guedira_seq_map/Allegro_test/groupA_mq20_region_bcfs"
mq_val=20


#### Executable ####

module load singularity/3.7.1
module load bcftools

echo
echo "Start time:"
date

regions_bed=$1
array_ind=$2


## Create output directory if array index is 0 or 1
if [[ $array_ind -eq 0 || $array_ind -eq 1 ]]; then
    mkdir -p "$out_dir"
    mkdir "${out_dir}/temp_files"
else
    sleep 10s
fi

## If the supplied array_ind is 0, we copy full regions .bed file to temp dir.
## Otherwise extract single region from .bed file
if [[ array_ind -eq 0 ]]; then
    label="all_regions"
    cp "$regions_bed" "${out_dir}/temp_files/${label}.bed"
else
    ## This ensures output files will have names in order, to avoid needing to
    ## sort after concattenating them together
    n=$(($(wc -l < "$regions_bed" | wc -c) - 1))
    prefix=$(printf "%0${n}d" $array_ind)

    region=$(head -n $array_ind "$regions_bed" | tail -n 1)
    suffix=$(echo "$region" | tr "\t" "_")
    label="${prefix}_${suffix}"
    echo "$region" > "${out_dir}/temp_files/${label}.bed"
fi

echo
echo "Input genomic region: ${label}"

## Perform the mpileup and calling
bcftools mpileup --fasta-ref "$ref_gen" \
                 --bam-list "$bams_list" \
                 --min-MQ "$mq_val" \
                 --regions-file "${out_dir}/temp_files/${label}.bed" \
                 --annotate FORMAT/DP,FORMAT/AD \
                 --output-type u |
   bcftools call --multiallelic-caller \
                 --variants-only \
                 --output-type b \
                 --output "${out_dir}/${label}.bcf"

echo
echo "End time:"
date
