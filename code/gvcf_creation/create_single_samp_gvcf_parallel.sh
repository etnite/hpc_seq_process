#!/bin/bash
set -e

## Create gVCF from single-sample BAM file
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## NOTE: Older versions of GATK cannot handle chromosomes as large as those found
## in wheat. Here, a Conda environment named "gatk" is used, as the GATK version
## installed on Ceres is not recent enough. This environment must contain the
## gatk4 package. To make it on Ceres, type:
##
##   module load miniconda
##   conda create --name gatk gatk4
##
## And go through the prompts.
##
## This script takes an input text file listing paths to .bam files, one per line
## It is designed to work with arrayer.sh, which will dispatch this script to run
## on one of the .bam files listed in the input text file. 
##
## WARNING: GATK's HaplotypeCaller is memory-hungry, and will easily overrun
## individual core memory limits if not reigned in. Memory limits should be
## set in the call to gatk below - for instance, setting:
##
##   gatk --java-options "-Xmx6g" HaplotypeCaller ...
##
## will limit memory usage to 6GB. Although HaplotypeCaller can output bgzipped
## gVCF files directly, this always seems to overrun memory in my testing. That
## is why this script uses the less efficient method of outputting a non-compressed
## VCF file, and then compressing and indexing it. Note that if this script is
## being run in parallel on many files simultaneously, this can lead to the
## temporary creation of files that take up a huge amount of space.
################################################################################


#### User-defined Constants ####

ref_gen="/project/genolabswheatphg/v1_refseq/whole_chroms/Triticum_aestivum.IWGSC.dna.toplevel.fa"
bam_list="/home/brian.ward/search_pattern_files/big_mq20_bams.txt"
out_dir="/lustre/project/guedira_seq_map/brian/US_excap/v1_mq20_gVCFs"


#### Executable ####

module load tabix
module load miniconda
source activate gatk4

mkdir -p "$out_dir"

## Read in array index (integer) - get corresponding .bam file name
## Generate output file name
arr_ind=$1
bam_file=$(head -n $arr_ind "$bam_list" | tail -n 1)
bam_base=$(basename "$bam_file")
samp="${bam_base%.*}"
out_file="${out_dir}/${samp}.g.vcf"

## Run haplotype caller
gatk --java-options "-Xmx3g" HaplotypeCaller \
     --reference "$ref_gen" \
     --input "$bam_file" \
     --read-index "${bam_file}.csi" \
     --emit-ref-confidence GVCF \
     --output "$out_file"

## Compress and index
bgzip --force "$out_file"
tabix --csi --preset vcf "${out_file}.gz" 
