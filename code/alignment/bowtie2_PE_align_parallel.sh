#!/bin/bash
shopt -s nullglob

## Paired end alignment using Bowtie2
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script aligns a single pair of mated fastq files, or one interleaved fastq file.
## It takes a single string as its only positional argument. Then it searches 
## for corresponding fastq files using the pattern fastq_dir/<pattern>*R[1/2]*[fastq/fq].gz
## for paired input, or fastq_dir/<pattern>*[fastq/fq].gz for interleaved.
## If the search pattern ends in "fastq.gz" or "fq.gz", then we assume the input
## is interleaved, since it's pointing to a single file.
##
## Everything up to the first underscore in the search string is used as the
## sample name.
##
## The output consists of a single sorted BAM file, with PCR duplicates marked.
##
## This script is intended to be used with arrayer.sh, to enable independent
## parallel runs on multiple samples simultaneously.
##
## NOTES: 
##
##   1) bowtie2 is run using the --sensitive-local option. This can be changed
##      manually in the call to bowtie2 below
##   2) The script will produce a .csi index of the output BAM file
##   3) This script requires two sorting steps in samtools, so can be a bit
##      slow...
##   4) Working directory inherited from parallelizing script - it is easiest
##      to define absolute paths
##   5) If the fastq file is from the sequence read archive (SRA), then its SRR
##      number is used as the barcode. Otherwise the barcode stored in the Illumina
##      read IDs is used.
##   6) Aligning 1 million paired reads (i.e. 2 million total reads) on Ceres
##      using one hyprethreaded core (2 threads) took 30 minutes.
##   7) Script assumes various read groups for a sample came from a single library
################################################################################


#### User-defined constants ####

## Reference genome fasta ("ref") must already be indexed using bowtie2-build
## and samtools index
fastq_dir="/project/genolabswheatphg/US_excap/filt_fastqs"
patterns_file="/home/brian.ward/search_pattern_files/ags_la_filt_fqs.txt"
out_dir="/project/genolabswheatphg/US_excap/v1_alignments"
ref="/project/genolabswheatphg/v1_refseq/whole_chroms/Triticum_aestivum.IWGSC.dna.toplevel.fa"

## Convert sample name in RG line of BAM header to uppercase? (TRUE/FALSE)
name2upper="TRUE"


#### Executable  ####

module load bowtie2
module load samtools

echo
echo "Start bowtie2_PE_align_parallel.sh"
echo "Start time:"
date

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
echo
echo "Input search pattern: ${patt}"

## Get sample name; convert sample name to uppercase if specified
samp=$(basename "$patt" | sed 's/_.*//')
if [[ "$name2upper" == "TRUE" ]]; then samp="${samp^^}"; fi

## Create output subdirectory from search pattern, if one is present
sub_dir=$(dirname "$patt")
if [[ "$subdir" != "." ]]; then mkdir -p "${out_dir}/${sub_dir}"; fi

## Generate output prefix. May need to remove up to two extensions, in case
## the search pattern ends with, e.g. ".fastq.gz"
out_pref="${out_dir}/${patt%.*}"
out_pref=$(echo "${out_pref%.*}" | sed 's/_interleaved//')


## If search pattern ends in fastq.gz or fq.gz, then we are assuming interleaved format
if [[ "$patt" == *fastq.gz ]] || [[ "$patt" == *fq.gz ]]; then
    fq="${fastq_dir}/${patt}"
else
    fq=$(echo "${fastq_dir}/${patt}"*R1*fastq.gz "${fastq_dir}/${patt}"*R1*fq.gz)
    fq2=$(echo "${fastq_dir}/${patt}"*R2*fastq.gz "${fastq_dir}/${patt}"*R2*fq.gz)
fi


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


## Run bowtie for either paired or interleaved formats
ref="${ref%.*}"
if [[ "$patt" == *fastq.gz ]] || [[ "$patt" == *fq.gz ]]; then
    bowtie2 -x "$ref" \
        --threads $SLURM_NTASKS \
        --rg-id "${samp}.${fcell}.${lane}.${bcode}" \
        --rg SM:"$samp" \
        --rg PL:ILLUMINA \
        --rg PU:"${fcell}.${lane}" \
        --rg LB:"${samp}.lib01" \
        --sensitive-local \
        --phred33 \
        --interleaved "$fq" |
        samtools sort -n -T "${out_pref}_sort1" -O SAM - |
        samtools fixmate -m -O SAM - - |
        samtools sort -T "${out_pref}_sort2" -O SAM - |
        samtools markdup - "${out_pref}.bam"
else
    bowtie2 -x "$ref" \
        --threads $SLURM_NTASKS \
        --rg-id "${samp}.${fcell}.${lane}.${bcode}" \
        --rg SM:"$samp" \
        --rg PL:ILLUMINA \
        --rg PU:"${fcell}.${lane}" \
        --rg LB:"${samp}.lib01" \
        --sensitive-local \
        --phred33 \
        -1 "$fq" \
        -2 "$fq2" |
        samtools sort -n -T "${out_pref}_sort1" -O SAM - |
        samtools fixmate -m -O SAM - - |
        samtools sort -T "${out_pref}_sort2" -O SAM - |
        samtools markdup - "${out_pref}.bam"
fi


## Index BAM file using .csi index format
samtools index -c "${out_pref}.bam"

echo
echo "End time:"
date
