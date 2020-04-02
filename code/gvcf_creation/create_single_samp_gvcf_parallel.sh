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
## This script assumes that BAM files are all located within a single directory,
## and that they labeled in the form: sample_name.bam
##
## WARNING: GATK's HaplotypeCaller is memory-hungry, and will easily overrun
## individual core memory limits if not reigned in. Memory limits should be
## set in the call to gatk below - for instance, setting:
##
##   gatk --java-options "-Xmx6g" HaplotypeCaller ...
##
## will limit memory usage to 6GB.
##
## In addition, this script is designed to work with arrayer.sh, to perform the 
## gVCF conversion on multiple samples simultaneously.
################################################################################


#### User-defined Constants ####

ref_gen="/project/genolabswheatphg/v1_refseq/whole_chroms/Triticum_aestivum.IWGSC.dna.toplevel.fa"
in_dir="/project/genolabswheatphg/alignments/ERSGGL_SRW_bw2_bams/SRW_merged_excap_GBS_wholechrom_bw2_bams_mq20_filt"
out_dir="/project/genolabswheatphg/gvcfs/SRW_single_samp_bw2_excap_GBS_mq20"


#### Executable ####

module load tabix
module load miniconda
source activate gatk

mkdir -p "${out_dir}"

## Read in array index (integer) - get corresponding .bam file name
## Generate output file name
arr_ind=$1
bam_file=$(ls -1 "${in_dir}"/*.bam | head -n $arr_ind | tail -n 1)
bam_base=$(basename "${bam_file}")
samp="${bam_base%.*}"
out_file="${out_dir}"/"${samp}".g.vcf

## Run haplotype caller
gatk --java-options "-Xmx3g" HaplotypeCaller \
     -R "${ref_gen}" \
     -I "${bam_file}" \
     --read-index "${bam_file}".csi \
     -ERC GVCF \
     -O "${out_file}"

## Compress and index
bgzip "${out_file}"
tabix --csi --preset vcf "${out_file}".gz 

source deactivate
