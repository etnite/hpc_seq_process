#!/bin/bash
set -e

## Filter VCF file - parallel
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script filters a VCF file variant-wise. It will also subset samples based
## on a user-defined list of sample names, if supplied.
##
## The script takes two positional, command-line arguments:
##   1) The path to a .bed file defining regions of VCF file to extract and filter.
##      This file must be sorted by chromosome and then start position, and should
##      not contain any overlapping regions.
##   2) An integer, defining which line of the .bed file to use. Supplying 0 will
##      use all lines of the .bed file. This integer will typically be supplied
##      by a parallel dispatch script, which will allow for filtering multiple
##      regions simultaneously.
##
## All other constants are set within the script, below this description.
##
## Whether the input is a VCF or BCF file, the output will always be a BCF file.
## If using all regions, the output file will be placed in the user-defined
## output directory as "all_regions.bcf". If a single region is filtered, then
## the output BCF file will have a name formatted as [region_integer]_[chrom]_[start]_[end].bcf
## The integer in the file name is formatted with leading zeroes to allow concatenation
## of multiple BCF files without a sorting step. For instance, an output BCF file name
## of region number 15 out of 1000 might be formatted as:
##
##   0015_2A_100000_200000.bcf
##
## NOTES:
##
##   The max_miss parameter works in an opposite manner of the VCFTools
##   implementation, which I find unintuitive. Here, if max_miss is set to 1, SNPs with 100%
##   missing data would theoretically be allowed, while setting max_miss to 0 will only
##   allow SNPs without any missing data
##
##   Samples in the output file will always be sorted.
##
##   Filtering to only retain biallelic variants, only retain SNPs, etc. must
##   be hand-coded into the first bcftools call in the piped command
##
##   The input VCF/BCF file must be indexed, and must include contig lines in
##   the header.
##
##   A regions .bed file is still required if filtering an entire VCF. This can
##   be created from the contig lines in the VCF header, e.g.:
##
##     bcftools view -h [input.vcf.gz] | grep "##contig" | sed -e 's/##contig=<ID=//' -e 's/,length=/\t0\t/' -e 's/[,>$]//'
##
##   For depth filtering, this script uses the DP value in the INFO column. This
##   value is the total depth of reads that passed the quality control parameters
##   defined during SNP calling, (e.g. minimum mapq value). FORMAT DP values, if
##   present, are for raw read counts per sample.
##
##   snpgap will remove SNPs within n bases of an indel (or overlapping)
##   indelgap will thin clusters of indels within n bases of each other, to only retain one
##
##   BCFTools is designed to work with .bcf files. If .vcf files are supplied, it
##   first converts them to .bcf format internally. Therefore, using a .bcf file
##   as input is faster, though BCF files may be larger than .vcf.gz files.
################################################################################


#### SLURM job control ####

#SBATCH --job-name="filt-vcf" #name of the job submitted
#SBATCH --partition=short #name of the queue you are submitting job to
  ##SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1  #Number of overall tasks - overrides tasks per node
  ##SBATCH --ntasks-per-node=22 #number of cores/tasks
#SBATCH --time=06:00:00 time allocated for this job hours:mins:seconds
#SBATCH --mail-user=bpward2@ncsu.edu #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### User-defined constants ####################################################

## Note that SNP depth and proportion of missing data are highly correlated

vcf_in="/lustre/project/genolabswheatphg/US_excap/v1_variants/brian_raw_v1_variants/v1_raw_variants.bcf"
out_dir="/lustre/project/genolabswheatphg/US_excap/v1_variants/brian_miss08_maf003_variants/v1_miss08_maf003_variants.bcf"
samp_file="none"
min_maf=0.03
max_miss=0.8
max_het=1
min_dp=0
max_dp=1e6
het2miss="true"
snpgap=3
indelgap=3


#### Executable ################################################################

## Using conda here to utilize bcftools plugins and plot-vcfstats
#module load bcftools
#module load miniconda
#source activate htslib

echo
echo "Start time:"
date

regions_bed=$1
array_ind=$2

## Get first letter of true/false het2miss option and
## Set some constants depending on het2miss
het2miss="${het2miss:0:1}"
if [[ "$het2miss" == [Tt] ]]; then
	max_het="NA"
	het_string="./."
fi

## Output filtering parameters and generate the samples file
if [[ $array_ind -eq 0 || $array_ind -eq 1 ]]; then
    mkdir -p "$out_dir"
    mkdir -p "${out_dir}/temp_files"

    ## Echo input parameters to output dir
    echo -e "Input VCF\t${vcf_in}" > "${out_dir}/filt_params.txt"
    echo -e "Regions .bed file\t${regions_bed}" >> "${out_dir}/filt_params.txt"
    echo -e "Sample subset list\t${samp_file}" >> "${out_dir}/filt_params.txt"
    echo -e "Minimum MAF\t${min_maf}" >> "${out_dir}/filt_params.txt"
    echo -e "Max missing proportion\t${max_miss}" >> "${out_dir}/filt_params.txt"
    echo -e "Max het. proportion\t${max_het}" >> "${out_dir}/filt_params.txt"
    echo -e "Min. average depth\t${min_dp}" >> "${out_dir}/filt_params.txt"
    echo -e "Max average depth\t${max_dp}" >> "${out_dir}/filt_params.txt"
    echo -e "Heterozygous calls converted to missing data?\t${het2miss}" >> "${out_dir}/filt_params.txt"
    echo -e "SNP overlap gap\t${snpgap}" >> "${out_dir}/filt_params.txt"
    echo -e "Indel overlap gap\t${indelgap}" >> "${out_dir}/filt_params.txt"

    ## If samp_file exists, use it to subset samples
    ## Otherwise retain all samples present in the input file
    if [[ -f "$samp_file" ]]; then
        cut -f 1 "$samp_file" | sort > "${out_dir}/temp_files/samp_file.txt"
    else
        bcftools query --list-samples "$vcf_in" | sort > "${out_dir}/temp_files/samp_file.txt"
    fi
else
    sleep 10s
fi

## If the supplied array_ind is 0, we copy full regions .bed file to temp dir.
## Otherwise extract single region from .bed file
if [[ array_ind -eq 0 ]]; then
    label="all_regions"
    cp "$regions_bed" "${out_dir}/temp_files/${label}.bed"
else
    ## This ensures output VCFs will have names in order, to avoid needing to
    ## sort after concattenating them together
    n=$(($(wc -l < "$regions_bed" | wc -c) - 1))
    prefix=$(printf "%0${n}d" $array_ind)

    region=$(head -n $array_ind "$regions_bed" | tail -n 1)
    suffix=$(echo "$region" | tr "\t" "_")
    label="${prefix}_${suffix}"
    echo "$region" > "${out_dir}/temp_files/${label}.bed"
fi



## Big copy-pasted if-else block for different actions depending het2miss
## Using het2miss will leave you with some SNPs that only have one allele, so some extra steps
## Are necessary to get rid of these
echo "Filtering VCF..."
echo

if [[ "$het2miss" == [Tt] ]]; then
    bcftools view "${vcf_in}" \
        --samples-file "${out_dir}/temp_files/samp_file.txt" \
        --min-alleles 2 \
        --max-alleles 2 \
        --output-type u |
    bcftools +setGT --output-type u - -- --target-gt q \
        --include 'GT="het"' \
        --new-gt "$het_string" |
    bcftools +fill-tags --output-type u -- -t MAF,F_MISSING |
    bcftools view - \
        --regions-file "${out_dir}/temp_files/${label}.bed" \
        --exclude "INFO/F_MISSING > ${max_miss} || INFO/MAF < ${min_maf} || INFO/DP < ${min_dp} || INFO/DP > ${max_dp} || (COUNT(GT=\"het\") / COUNT(GT!~\"\.\")) > ${max_het}" \
        --output-type u |
    bcftools filter --SnpGap $snpgap \
        --IndelGap $indelgap \
        --output-type b \
        --output "${out_dir}/${label}.bcf"
else
    bcftools view "${vcf_in}" \
        --samples-file "${out_dir}/temp_files/samp_file.txt" \
        --min-alleles 2 \
        --max-alleles 2 \
        --output-type u |
    bcftools +fill-tags --output-type u -- -t MAF,F_MISSING |
    bcftools view - \
        --regions-file "${out_dir}/temp_files/${label}.bed" \
        --exclude "INFO/F_MISSING > ${max_miss} || INFO/MAF < ${min_maf} || INFO/DP < ${min_dp} || INFO/DP > ${max_dp}" \
        --output-type u |
    bcftools filter --SnpGap $snpgap \
        --IndelGap $indelgap \
        --output-type b \
        --output "${out_dir}/${label}.bcf"
fi

## If a whole VCF/BCF file was filtered, generate some summary stats and plots
if [[ $array_ind -eq 0 ]]; then
    bcftools index -c "${out_dir}/${label}.bcf"

    ## Generate summary stats
    bcftools stats "$vcf_out" > "${out_dir}/stats.txt"

    ## Create plots from summary stats
    mkdir "${out_dir}/plots"
    plot-vcfstats --prefix "${out_dir}/plots" --no-PDF "${out_dir}/stats.txt"
fi

## Generate summary stats using TASSEL
#echo "Generating summary statistics with TASSEL..."
#$TASSEL_PL -vcf "${vcf_out}" \
#           -GenotypeSummaryPlugin \
#           -endPlugin \
#           -export "${out_dir}"/summary

rm -rf "${out_dir}/temp_files"

echo
echo "End time:"
date
