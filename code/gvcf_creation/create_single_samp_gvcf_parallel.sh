#!/bin/bash
set -e

## Create gVCF from single-sample BAM file
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
## This script assumes that BAM files are all located within a single directory,
## and that they labeled in the form: sample_name.bam
##
## In addition, this script is designed to work with arrayer.sh, to perform the 
## gVCF conversion on multiple samples simultaneously.
################################################################################


#### User-defined Constants ####

ref_gen="/project/genolabswheatphg/v1_refseq/whole_chroms/Triticum_aestivum.IWGSC.dna.toplevel.fa"
in_dir="/project/genolabswheatphg/alignments/SRW_filt_bams"
out_dir="/project/genolabswheatphg/gvcfs/SRW_gvcfs"


#### Executable ####

module load tabix
module load miniconda
source activate gatk


## Read in array index (integer) - get corresponding .bam file name
## Generate output file name
arr_ind=$1
bam_file=$(ls -1 "${in_dir}"/*.bam | head -n $arr_ind | tail -n 1)
bam_base=$(basename "${bam_file}")
samp="${bam_base%.*}"
out_file="${out_dir}"/"${samp}".g.vcf

## Ran haplotype caller
gatk HaplotypeCaller \
    -R "${ref_gen}" \
    -I "${bam_file}" \
    --read-index "${bam_file}".csi \
    -ERC GVCF \
    -O "${out_file}"

## Compress and index
bgzip "${out_file}"
tabix --csi --preset vcf "${out_file}".gz 

source deactivate

