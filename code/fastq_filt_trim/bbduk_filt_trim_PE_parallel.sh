#!/bin/bash
shopt -s nullglob

## Paired-end FASTQ quality filtering and adapter trimming using BBDuk
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This version of the script runs BBDuk on a single pair of mated fastq files.
## It takes a single string (usually a sample name) as its only positional
## argument. Then it searches for corresponding fastq files using the pattern:
## input_dir/samplename_R[12].fastq.gz
##
## This script is intended to be used with parallelize.sh, to enable independent
## parallel runs on multiple samples simultaneously.
##
## NOTES: 
##   1) bbduk parameters are hard coded in the bbduk.sh call below. Note -Xmx
##      argument controlling memory used by bbduk
##   2) Working directory inherited from parallelizing script - it is easiest
##      to define absolute paths
##   3) Processing ~50 million paired end reads with one core took around
##      3 mins. on Ceres cluster
################################################################################


#### User-defined constants ####

in_dir="/project/genolabswheatphg/merged_fastqs/v1_hapmap"
adapt_fasta="../../bbduk_38.84_Illumina_adapters.fa"
patterns_file="/home/brian.ward/repos/wheat_phg/sample_lists/v1_hapmap_bioproj/sample_names.txt"
out_dir="/project/genolabswheatphg/filt_fastqs/v1_hapmap"


#### Executable  ####

module load bbtools

echo
echo "Start bbduk_filt_trim_paired_parallel.sh"
echo "Start time:"
date

mkdir -p "${out_dir}"
array_ind=$1

## Get length of adapters
#ad_len=$(head -n 2 "${adapt_fasta}" | tail -n -1 | wc -c)

## Get search pattern
patt=$(head -n "${array_ind}" "${patterns_file}" | tail -n 1)

## Get sub-directory from search pattern, if present
sub_dir=$(dirname "$patt")
if [[ "$subdir" != "." ]]; then mkdir -p "${out_dir}/${sub_dir}"; fi

## If search pattern ends in fastq.gz or fq.gz, then we are assuming interleaved format
if [[ "$patt" == *fastq.gz ]] || [[ "$patt" == *fq.gz ]]; then
    fq="${in_dir}/${patt}"

    ## Run BBDuk
    ## Setting maq will remove entire reads based on their average qual
    ## maq=10 average error is 10%
    ## maq=13 ~5%
    ## maq=15 ~3%
    ## maq=20 1%
    bbduk.sh -Xmx2800m in="${fq}" \
        out="${out_dir}/${patt}" \
        ref="$adapt_fasta" ktrim=r k=23 mink=11 hdist=1 ftm=5 maq=13 minlen=75 tpe tbo

else
    fq=$(echo "${in_dir}"/"${patt}"*R1*fastq.gz "${in_dir}"/"${patt}"*R1*fq.gz)
    fq2=$(echo "${in_dir}"/"${patt}"*R2*fastq.gz "${in_dir}"/"${patt}"*R2*fq.gz)

    ## Run BBDuk
    bbduk.sh -Xmx2800m in1="$fq" in2="$fq2" \
        out1="${out_dir}/${patt}_R1.fastq.gz" out2="${out_dir}/${patt}_R2.fastq.gz" \
        ref="$adapt_fasta" ktrim=r k=23 mink=11 hdist=1 ftm=5 maq=13 minlen=75 tpe tbo
fi

echo
echo "End time:"
date
