#!/bin/bash

#### WARNING: Under Development ####

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
in_dir="/autofs/bioinformatics-ward/MRAseq_v2_test_set/ubams"
mark_dir="/autofs/bioinformatics-ward/MRAseq_v2_test_set/marked_ubams"
out_dir="/autofs/bioinformatics-ward/MRAseq_v2_test_set/bams"
temp_dir="/scratch"
bams_file="/home/ward.1660/pattern_files/ubams_list.txt"
nthreads=4

ref="/autofs/bioinformatics-ward/CSv1_ref_chromparts/161010_Chinese_Spring_v1.0_pseudomolecules_parts"


#### Executable ####

#module load bowtie2
#module load bwa
#module load samtools
#source activate gatk4

echo
echo "ubam_to_bam.sh"
echo "Start time:"
date

array_ind=$1
mkdir -p "${out_dir}"
mkdir -p "$mark_dir"

## Get name of input BAM file
in_bam=$(head -n "${array_ind}" "${bams_file}" | tail -n 1)
#in_bam="TRIBUTE_sub10K_L001.bam"
mark_bam=$(echo "$in_bam" | sed 's/.bam/markadapters.bam/')
met_file=$(echo "$in_bam" | sed 's/.bam/markadapters_metrics.txt/')
out_bam=$(echo "$in_bam" | sed 's/.bam/_aligned.bam/')

## Mark Illumina adapters in uBAM
gatk MarkIlluminaAdapters --java-options "-Xmx12G" \
	-I "${in_dir}/${in_bam}" \
	-O "${mark_dir}/${mark_bam}" \
	-M "${mark_dir}/${met_file}" \
	--TMP_DIR "$temp_dir"

## Pipe to run SamtoFastq, followed by alignment, followed by MergeBamAlignment
gatk SamToFastq --java-options "-Xmx12G" \
    -I "$mark_dir"/"$mark_bam" \
    --FASTQ /dev/stdout \
    --CLIPPING_ATTRIBUTE XT \
    --CLIPPING_ACTION 2 \
    --INTERLEAVE false \
    --INCLUDE_NON_PF_READS true \
    --TMP_DIR "$temp_dir" |

#bwa mem -M -t $nthreads -p "${ref}" /dev/stdin |
bowtie2 -x "$ref" --threads $nthreads --sensitive-local --phred33 -U - > "${out_dir}/${out_bam}"


exit 0;

gatk MergeBamAlignment --java-options "-Xmx12G" \
    ALIGNED_BAM=/dev/stdin \
    UNMAPPED_BAM=$"{mark_dir}/{mark_bam}" \
    OUTPUT=$"{out_dir}/{out_bam}" \
    R="${ref}.fasta" \
    CREATE_INDEX=true \
    ADD_MATE_CIGAR=true \
    CLIP_ADAPTERS=false \
    CLIP_OVERLAPPING_READS=true \
    INCLUDE_SECONDARY_ALIGNMENTS=true \
    MAX_INSERTIONS_OR_DELETIONS=-1 \
    PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
    ATTRIBUTES_TO_RETAIN=XS \
    TMP_DIR="$temp_dir"    
