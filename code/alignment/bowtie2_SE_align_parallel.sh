#!/bin/bash

## Single end alignment using Bowtie2
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script aligns the reads in a single fastq file.
## It takes a single string (usually a sample name) as its only positional
## argument. Then it searches for the corresponding fastq file using the pattern:
## fastq_dir/samplename.fastq.gz
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
################################################################################


#### User-defined constants ####

## Reference genome fasta ("ref") must already be indexed using bowtie2-build
## and samtools index
#fastq_dir="/project/genolabswheatphg/filt_fastqs/SRW_excap"
#samps_file="/home/brian.ward/repos/wheat_phg/sample_lists/SRW_reform_samples.txt"
#out_dir="/project/genolabswheatphg/alignments/SRW_wholechrom_bw2_bams"
#ref="/project/genolabswheatphg/v1_refseq/whole_chroms/Triticum_aestivum.IWGSC.dna.toplevel.fa"

#fastq_dir="/project/genolabswheatphg/filt_fastqs/KS_HRW_excap"
#samps_file="/home/brian.ward/repos/wheat_phg/sample_lists/ZENDA.txt"
#out_dir="/project/genolabswheatphg/alignments/KS_HRW_wholechrom_bw2_bams"
#ref="/project/genolabswheatphg/v1_refseq/whole_chroms/Triticum_aestivum.IWGSC.dna.toplevel.fa"

#fastq_dir="/project/genolabswheatphg/filt_fastqs/wheatCAP_parents"
#samps_file="/home/brian.ward/repos/wheat_phg/sample_lists/wheatCAP_samples_reform.txt"  
#out_dir="/project/genolabswheatphg/alignments/wheatCAP_wholechrom_bw2_bams" 
#ref="/project/genolabswheatphg/v1_refseq/whole_chroms/Triticum_aestivum.IWGSC.dna.toplevel.fa" 

#fastq_dir="/project/genolabswheatphg/merged_fastqs/SRW_GBS"
#samps_file="/home/brian.ward/repos/wheat_phg/sample_lists/SRW_GBS_samps.txt"
#out_dir="/project/genolabswheatphg/alignments/SRW_GBS_wholechrom_bw2_bams"
#ref="/project/genolabswheatphg/v1_refseq/whole_chroms/Triticum_aestivum.IWGSC.dna.toplevel.fa"

fastq_dir="/project/genolabswheatphg/merged_fastqs/SRW_GBS"
samps_file="/home/brian.ward/repos/wheat_phg/sample_lists/SS-MPV57.txt"
out_dir="/project/genolabswheatphg/alignments/ERSGGL_SRW_bw2_bams/SRW_GBS_wholechrom_bw2_bams"
ref="/project/genolabswheatphg/v1_refseq/whole_chroms/Triticum_aestivum.IWGSC.dna.toplevel.fa"


#### Executable  ####

module load bowtie2
module load samtools

echo
echo "Start bowtie2_align_parallel.sh"
echo "Start time:"
date

mkdir -p "${out_dir}"
array_ind=$1

## Get sample name
samp=$(head -n "${array_ind}" "${samps_file}" | tail -n 1)

## Set alias for fastq file
fq=$(echo "${fastq_dir}"/"${samp}".fastq.gz)

## For some reason, the barcode indexes in FASTQ files can contain some N
## values for the first few reads. I don't know the significance of this. 
## Let's just grab line 10,001:
## NOTE: In the case of concatenated fastq files, flowcells, lanes,
## barcodes, etc. may all vary. Therefore these sections are disabled
#one_line=$(zcat $fq1 | head -n 10001 | tail -n -1)

## This sets the first four fields (instrument:run_id:flowcell:lane) of the 
## first line of the fastq file as a variable named "id"
#id=$(echo "$one_line" | cut -f 1-4 -d":" | sed 's/@//' | sed 's/:/_/g')

## Now get the barcode from the 10th field of the line
## NO
#bar=$(echo "$one_line" | cut -f 10 -d":" | sed 's/+/-/')


ref="${ref%.*}"
bowtie2 -x "${ref}" \
        --threads $SLURM_NTASKS \
        --rg-id "${samp}" \
        --rg SM:"${samp}" \
        --rg PL:ILLUMINA \
        --sensitive-local \
        --phred33 \
        -U "${fq}" |
        samtools sort -T "${out_dir}"/"${samp}"sort1 -O BAM -o "${out_dir}"/"${samp}".bam

## Index BAM file using .csi index format
samtools index -c "${out_dir}"/"${samp}".bam

echo
echo "End time:"
date
