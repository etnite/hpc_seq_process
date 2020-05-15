#!/bin/bash

## Convert a uBAM file to a BAM
##
## This script is designed to take an unmapped BAM file, and then:
##   1) Mark the positions of Illumina adapters within reads
##   2) align reads using bowtie2 (note that GATK best practices uses BWA MEM
##      for alignment)
##   3) Merge read info from the unaligned BAM into the aligned BAM
################################################################################


#### User-defined variables ####

## Input fastq directory, output directory for .bam files, and file
## of sample names, one listed per line
in_dir="/project/genolabswheatphg/merged_fastqs/SRW_excap/fq_demux_test"
mark_dir="/project/genolabswheatphg/merged_fastqs/SRW_excap/fq_demux_test/marked_adapters"
#out_dir="/project/genolabswheatphg/merged_fastqs/SRW_excap/fq_demux_test"
#bams_file="/home/brian.ward/repos/wheat_phg/sample_lists/v1_hapmap_bioproj/sample_names.txt"


#### Executable ####

module load bowtie2
module load samtools
source activate gatk4

echo
echo "ubam_to_bam.sh"
echo "Start time:"
date

array_ind=$1
mkdir -p "${out_dir}"
mkdir -p "$mark_dir"

## Get name of input BAM file
#in_bam=$(head -n "${array_ind}" "${bams_file}" | tail -n 1)
in_bam="TRIBUTE_sub10K_L001.bam"
mark_bam=$(echo "$in_bam" | sed 's/.bam/markadapters.bam/')
met_file=$(echo "$in_bam" | sed 's/.bam/markadapters_metrics.txt/')

## Mark Illumina adapters in uBAM
gatk MarkIlluminaAdapters \
	-I "$in_dir"/"$in_bam" \
	-O "$mark_dir"/"$mark_bam" \
	-M "$mark_dir"/"$met_file" \
	TMP_DIR="$mark_dir"

source deactivate
echo
echo "End time:"
date

exit 0;

## Pipe to run SamtoFastq, followed by alignment, followed by MergeBamAlignment

gatk -Xmx8G -jar SamToFastq \
I=6483_snippet_markilluminaadapters.bam \
FASTQ=/dev/stdout \
CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
TMPDIR=/path/shlee | \
/path/bwa mem -M -t 7 -p /path/Homosapiens_assembly19.fasta /dev/stdin | \


java -Xmx16G -jar /path/picard.jar MergeBamAlignment \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM=6383_snippet_revertsam.bam \
OUTPUT=6483_snippet_piped.bam \
R=/path/Homo_sapiens_assembly19.fasta CREATE_INDEX=true ADD_MATE_CIGAR=true \
CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
TMP_DIR=/path/shlee