#!/bin/bash
#set -e

## Rename samples, annotate variants, split SNPs and indels
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## INPUTS:
##
##   1. VCF file
##   2. reference .fasta file (must be indexed with samtools)
##   3. .gff3 annotations file
##   4. (optional) two-column tab-delimited file, with current sample names in
##      first column, and corresponding desired sample names in second column. 
##   5. "true" or "false" supplied for "split_file" argument, which controls
##      whether or not indels and SNPs are separated into different output files
##
## OUTPUTS:
##
##   1. Input VCF file, with variant consequences added, and (possibly) with
##      sample names updated.
##   2. New VCF file containing just SNPs
##   3. New VCF file containing just indels
################################################################################


#### SLURM job control #### 

#SBATCH --job-name="annote-vcf" #name of the job submitted
#SBATCH --partition=short #name of the queue you are submitting job to
  ##SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1  #Number of overall tasks - overrides tasks per node
  ##SBATCH --ntasks-per-node=22 #number of cores/tasks
#SBATCH --time=08:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=bpward2@ncsu.edu #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### User-defined constants ####

vcf_in="/lustre/project/genolabswheatphg/US_excap/v1_variants/brian_raw_v1_variants/v1_raw_variants.bcf"
ref="/lustre/project/genolabswheatphg/v1_refseq/whole_chroms/Triticum_aestivum.IWGSC.dna.toplevel.fa"
gff="/lustre/project/genolabswheatphg/v1_refseq/whole_chroms/Triticum_aestivum.IWGSC.44.gff3.gz"
lookup_file="none"
split_file="false"


#### Executable ####

module load bcftools

echo
echo "Start time:"
date

## Isolate the file format and extension
out_dir=$(dirname "$vcf_in")
ext="${vcf_in##*.}"
base="${vcf_in%.*}"
base="${vcf_in%.*}"
if [[ "$ext" == "gz" ]]; then
    fmt="z"
    ext="vcf.gz"
elif [[ "$ext" == "bcf" ]]; then
    fmt="b"
else
    echo
    echo "Input file should be either in BCF (.bcf) or gzipped VCF (.vcf.gz) format"
    exit 1;
fi


## If sample name lookup file is supplied, then replace sample names
## and call consequences. Otherwise, skip sample renaming step
if [[ -f "$lookup_file" ]]; then
	bcftools reheader "$vcf_in" --samples "$lookup_file" |
        bcftools csq -f "$ref" -g "$gff" -p a -O "$fmt" -o "${base}_csq.${ext}"
else
    bcftools csq "$vcf_in" -f "$ref" -g "$gff" -p a -O "$fmt" -o "${base}_csq.${ext}"
fi

bcftools index -c "${base}_csq.${ext}"

## Create SNP-only and indel-only files if the user specified
split_file="${split_file:0:1}"
split_file="${split_file^^}"
if [[ "$split_file" == "T" ]]; then
    bcftools view "${base}_csq.${ext}" --type snps -O "$fmt" -o "${base}_csq_snps_only.${ext}"
    bcftools index -c "${base}_csq_snps_only.${ext}"
    bcftools view "${base}_csq.${ext}" --type indels -O "$fmt" -o "${base}_csq_indels_only.${ext}"
    bcftools index -c "${base}_csq_indels_only.${ext}"
else
    echo
    echo "User selected not to split SNPs and indels into separate files"
fi

echo
echo "End time:"
date
