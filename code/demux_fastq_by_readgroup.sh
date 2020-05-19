#!/bin/bash
shopt -s nullglob

## Demultiplex paired fastq files by read group
##
## This script is designed to take a pair of fastq files, and demultiplex them
## by their read groups, outputting either interleaved fastq files, or else
## interleaved unaligned BAM (uBAM) files.
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
## This script takes as input either a single pair of mated fastq files, or one 
## interleaved fastq file. It takes a single string as its only positional argument. 
## Then it searches for corresponding fastq files using the pattern: 
## fastq_dir/<pattern>*R[1/2]*[fastq/fq].gz for paired input, or 
## fastq_dir/<pattern>*[fastq/fq].gz for interleaved.
## If the search pattern ends in "fastq.gz" or "fq.gz", then we assume the input
## is interleaved, since it's pointing to a single file.
##
## A search pattern can contain a subdirectory prefix. In this case, the
## subdirectory will be created in the output directory.
##
## BBDuk's demuxbyname.sh is used to split apart the fastq files and write out
## a single interleaved fastq for each read group. Finally, Picard FastqToSam
## converts each fastq to a uBAM, adding read group metadata.
################################################################################



#### User-defined variables ####

## Input fastq directory, output disrectory for .fastq or .bam files, file
## of file name patterns, one listed per line, and desired output format, either
## "fq" or "bam"
fastq_dir="/project/genolabswheatphg/merged_fastqs/SRW_excap"
out_dir="/project/genolabswheatphg/merged_fastqs/SRW_excap/fq_demux_test"
patterns_file="/home/brian.ward/test_samplist2.txt"
out_fmt="bam"

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


## Sanity check on output format
if [[ "$out_fmt" != "fq" ]] && [[ "$out_fmt" != "bam" ]]; then
    echo
    echo "Error - please set output format (out_fmt) to either 'fq' or 'bam'"
    echo "It is currently set to: ${out_fmt}"
    exit 1;
fi

mkdir -p "$out_dir"
array_ind=$1

## Get search pattern string and sample name; convert sample name to uppercase
#patt="TRIBUTE_ileaved_sub10K.fastq.gz"
patt=$(head -n "$array_ind" "$patterns_file" | tail -n 1)
samp=$(basename "$patt" | sed 's/_.*//')
upsamp="${samp^^}"
sub_dir=$(dirname "$patt")
if [[ "$subdir" != "." ]]; then
	mkdir -p "${out_dir}/${sub_dir}"
	upsamp="${sub_dir}/${upsamp}"
fi


## If search pattern ends in fastq.gz or fq.gz, then we are assuming interleaved format
if [[ "$patt" == *fastq.gz ]] || [[ "$patt" == *fq.gz ]]; then
    fq="${fastq_dir}/${patt}"
else
    fq=$(echo "${fastq_dir}/${patt}"*R1*fastq.gz "${fastq_dir}/${patt}"*R1*fq.gz)
    fq2=$(echo "${fastq_dir}/${patt}"*R2*fastq.gz "${fastq_dir}/${patt}"*R1*fq.gz)
fi


## Read every nth line of the first read fastq file
## Need if statement to check for SRR files
zcat "$fq" | sed -n "1~${read_samp}p" | cut -d ":" -f 3,4 | sort -u > "${out_dir}/${upsamp}_fcell_lane.txt"


## Demultiplex the fastq based on the unique flowcell lane combination
if [[ "$patt" == *fastq.gz ]] || [[ "$patt" == *fq.gz ]]; then
    demuxbyname.sh in="$fq" \
    	substringmode \
    	out="${out_dir}/${upsamp}"_%_interleaved.fastq.gz \
    	names="${out_dir}/${upsamp}_fcell_lane.txt"
else
    demuxbyname.sh in1="$fq" in2="$fq2" \
    	substringmode \
    	out="${out_dir}/${upsamp}"_%_interleaved.fastq.gz \
    	names="${out_dir}/${upsamp}_fcell_lane.txt"
fi


## Loop through new demuxed, interleaved fastq files
out_fqs=( "${out_dir}/${upsamp}"*interleaved.fastq.gz )
for i in "${out_fqs[@]}"; do

	if [[ "$out_fmt" == "fq" ]]; then

		## The interleaved fastq has a colon between flowcell and lane
		## Replace with underscore
		new_name=$(echo "$i" | sed 's/:/_/g')
		mv "$i" "$new_name"

	else

	    ## Get 10,001th line. Sometimes barcodes can contain Ns at beginning of file
        ## Test for whether dealing with SRA files
        id_line=$(zcat "$i" | head -n 10001 | tail -n 1)
        fcell=$(echo "$id_line" | cut -d ":" -f 3)
        lane=$(echo "$id_line" | cut -d ":" -f 4)
        if [[ "$id_line" == @SRR* ]]; then
            bcode=$(echo "$id_line" | cut -d "." -f 1)
        else
            bcode=$(echo "$id_line" | cut -d ":" -f 10)
        fi

	    ## Create output .bam file name
	    out_bam=$(echo "$i" | sed 's/:/_/g' |  sed 's/.fastq.gz/.bam/')

	    ## Convert fastq to bam, adding metadata
	    gatk FastqToSam \
	        -F1 "$i" \
	        -O "$out_bam" \
	        -RG "${fcell}_${lane}" \
	        -SM $(basename "$upsamp") \
	        -PL ILLUMINA \
	        -PU "${fcell}_${lane}.${bcode}" \
	        -LB "${upsamp}_lib"

	    rm "$i"
	fi
done

rm "${out_dir}/${upsamp}_fcell_lane.txt"

source deactivate
echo
echo "End time:"
date
