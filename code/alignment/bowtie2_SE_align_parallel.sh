#!/bin/bash

## Single end alignment using Bowtie2
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script aligns the reads in a single fastq file.
## It takes a single string (a fastq file name) as its only positional argument. 
## Then it searches for the corresponding fastq file using the pattern:
## <fastq_dir>/<pattern>
##
## Everything up to the first underscore in the search string is used as the
## sample name.
##
## The output consists of a single sorted BAM file.
##
## This script is intended to be used with arrayer.sh, to enable independent
## parallel runs on multiple samples simultaneously.
##
## NOTES: 
##
##   1) bowtie2 is run using the --sensitive-local option. This can be changed
##      manually in the call to bowtie2 below
##   2) The script will produce a .csi index of the output BAM file
##   3) This script requires a sorting step in samtools, so can be a bit
##      slow...
##   4) Working directory inherited from parallelizing script - it is easiest
##      to define absolute paths
##   5) If the fastq file is from the sequence read archive (SRA), then its SRR
##      number is used as the barcode. Otherwise the barcode stored in the Illumina
##      read IDs is used.
##   6) Script assumes various read groups for a sample came from a single library
################################################################################


#### User-defined constants ####

## Reference genome fasta ("ref") must already be indexed using bowtie2-build
## and samtools index
fastq_dir="/project/guedira_seq_map/Allegro_test/filt_fastq"
patterns_file="/home/brian.ward/samp_and_file_lists/Allegro_filt_fqs_list.txt"
out_dir="/project/guedira_seq_map/Allegro_test/bams"
ref="/project/guedira_seq_map/ref_genomes/v1_refseq_w_KIMs/CSv1_refseq_w_KIMs.fa"

## Convert sample name to uppercase? (TRUE/FALSE) 
name2upper="TRUE"

## Number of threads - set to $SLURM_NTASKS to use number of available cores defined by scheduler
nthreads=$SLURM_NTASKS


#### Executable  ####

module load singularity/3.7.1
module load bowtie2
module load samtools

echo
echo "Start bowtie2_SE_align_parallel.sh"
echo "Start time:"
date

echo
echo "Number of cores: ${nthreads}"

## Input sanity check
name2upper="${name2upper^^}"
if [[ "$name2upper" != "TRUE"  ]] && [[ "$name2upper" != "FALSE" ]]; then
    echo
    echo "Error - please set name2upper to either 'TRUE' or 'FALSE'"
    echo "It is currently set to: ${name2upper}"
    exit 1;
fi

mkdir -p "$out_dir"
array_ind=$1

## Get search pattern string
patt=$(head -n "$array_ind" "$patterns_file" | tail -n 1)

## Get sample name; convert sample name to uppercase if specified
samp=$(basename "$patt" | sed 's/_.*//')
if [[ "$name2upper" == "TRUE" ]]; then samp="${samp^^}"; fi

## Create output subdirectory from search pattern, if one is present
sub_dir=$(dirname "$patt")
if [[ "$subdir" != "." ]]; then mkdir -p "${out_dir}/${sub_dir}"; fi

## Generate output prefix. May need to remove up to two extensions, in case
## the search pattern ends with, e.g. ".fastq.gz"
out_pref="${out_dir}/${patt%.*}"
out_pref="${out_pref%.*}"

## Set the full path to the fastq file
fq="${fastq_dir}/${patt}"


## Get RG tag info from fastq read IDs
## For some reason, the barcode indexes in FASTQ files can contain some N
## values for the first few reads. I don't know the significance of this. 
## Let's just grab line 50,001:
## NOTE: Concatenated fastqs will contain multiple read groups
id_line=$(zcat "$fq" | head -n 50001 | tail -n 1)
fcell=$(echo "$id_line" | cut -d ":" -f 3)
lane=$(echo "$id_line" | cut -d ":" -f 4)
if [[ "$id_line" == @SRR* ]]; then
    bcode=$(echo "$id_line" | cut -d "." -f 1 | sed 's/^@//')
else
    bcode=$(echo "$id_line" | cut -d ":" -f 10)
fi


## Run bowtie2
ref="${ref%.*}"
bowtie2 -x "$ref" \
        --threads $nthreads \
        --rg-id "${samp}.${fcell}.${lane}.${bcode}" \
        --rg SM:"$samp" \
        --rg PL:ILLUMINA \
        --rg PU:"${fcell}.${lane}" \
        --rg LB:"${samp}.lib01" \
        --sensitive-local \
        --phred33 \
        -U "$fq" |
        samtools sort -T "${out_pref}_sort1" -O BAM -o "${out_pref}.bam"

## Index BAM file using .csi index format
samtools index -c "${out_pref}.bam"

echo
echo "End time:"
date
