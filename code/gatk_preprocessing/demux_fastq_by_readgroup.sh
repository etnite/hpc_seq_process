#!/bin/bash

## Demultiplex paired fastq files by read group and convert to uBAMs
##
## This script is designed to take a pair of fastq files, and demultiplex them
## by their read groups, finally converting them to interleaved unaligned BAM
## (uBAM) files.
##
## uBAM files are  the preferred storage medium for unaligned read data in the
## GATK best practices pipeline. Specifically, there should be one uBAM per
## read group. What is a read group? Essentially, it is one run of a sequencing
## instrument. For instance, if a single sample is sequenced across multiple
## Illumina lanes, then these are separate read groups, despite being prepared
## from the same sample and library.
##
## This script is specifically for use when we have fastq files containing more
## than one read group, but we don't know what these read groups are ahead of
## time. This can happen with data downloaded from the NCBI Sequence Read Archive,
## or the European Nucleotide Archive. It can also happen when someone gives
## us single-sample fastq files that consist of multiple read-group fastq files
## concatenated together.
##
## Unfortunately in this case we must use a brute-force method to find the read
## groups, by reading every nth line (i.e. the read ID lines) of one of the paired 
## files, parsing flowcell, lane, and barcode info, and then retaining only unique
## values. The script WILL ONLY WORK WITH 4-LINE fastqs. If you have fastqs that
## split read strings across multiple lines, that will lead to big problems
##
## The script assumes that a sample has a single pair of fastq files, with
## names in the form <SAMPLE>*R1.fastq.gz and <SAMPLE>*R2.fastq.gz respectively.
##
## BBDuk's demuxbyname.sh is used to split apart the fastq files and write out
## a single interleaved fastq for each read group. Finally, Picard FastqToSam
## converts each fastq to a uBAM, adding read group metadata.
################################################################################



#### User-defined variables ####

## Input fastq directory, output directory for .bam files, and file
## of sample names, one listed per line
fastq_dir="/project/genolabswheatphg/merged_fastqs/SRW_excap"
out_dir="/project/genolabswheatphg/merged_fastqs/SRW_excap/fq_demux_test"
#samp_file="/home/brian.ward/repos/wheat_phg/sample_lists/v1_hapmap_bioproj/sample_names.txt"

## How often to sample read names from each fastq file
## For instance, setting to 10000 will sample every 10,000th read to
## try and find unique flowcell/lane combinations
read_samp=250


#### Executable ####

module load bbtools
source activate gatk4

echo
echo "demult_fastq_by_readgroup.sh"
echo "Start time:"
date

array_ind=$1
mkdir -p "${out_dir}"

## Get sample name
#samp=$(head -n "${array_ind}" "${samp_file}" | tail -n 1)
samp="TRIBUTE_sub10K"

## Set forward and reverse read fastq files
fq1=$(echo "${fastq_dir}"/"${samp}"*R1*fastq.gz)
fq2=$(echo "${fastq_dir}"/"${samp}"*R2*fastq.gz)

## Convert sample name to uppercase and underscores to dashes
## Personal pref - I like samples to have dashes to make pattern recognition easier
upsamp="${samp^^}"
upsamp=$(echo "${upsamp}" | sed 's/_/-/g')

## Read every nth line of the first read fastq file (typically 10,000)
## Need if statement to check for SRR files
zcat "$fq1" | sed -n "1~${read_samp}p" | cut -d ":" -f 3,4 | sort -u > "$out_dir"/fcell_lane.txt

## Demultiplex the fastq based on the unique flowcell lane combination
demuxbyname.sh in1="$fq1" in2="$fq2" \
    substringmode \
    out="$out_dir"/"$upsamp"_%.fastq.gz \
    names="$out_dir"/fcell_lane.txt

## Loop through new demuxed, interleaved fastq files
out_fqs=( "$out_dir"/"$upsamp"*.fastq.gz )
for i in "${out_fqs[@]}"; do

    ## Get 10,001th line. Sometimes barcodes can contain Ns at beginning of file
    id_line=$(zcat "$i" | head -n 10001 | tail -n 1)
    fcell=$(echo $id_line | cut -d ":" -f 3)
    lane=$(echo $id_line | cut -d ":" -f 4)
    bcode=$(echo $id_line | cut -d ":" -f 10)

    ## The interleaved fastq has a colon between flowcell and lane
    ## This is replaced with an underscore in the .bam file
    out_bam=$(echo "$i" | sed 's/:/_/g' |  sed 's/.fastq.gz/.bam/')

    gatk FastqToSam \
        -F1 "$i" \
        -O "$out_bam" \
        -RG "$fcell"_"$lane" \
        -SM "$upsamp" \
        -PL ILLUMINA \
        -PU "$fcell"_"$lane"."$bcode" \
        -LB "$upsamp"_lib

    rm "$i"
done

source deactivate
echo
echo "End time:"
date
