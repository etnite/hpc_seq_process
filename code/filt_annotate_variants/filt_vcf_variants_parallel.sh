#!/bin/bash
set -e


## Filter VCF file variant-wise and sample-wise - single threaded
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script filters a VCF file, both variant-wise and sample-wise. It does this
## in the following order:
##   1. Filter variants according to thresholds and retain only user-specified samples
##   2. Optionally filter samples using thresholds
##
## NOTE: The max_miss parameter works in an opposite manner of the VCFTools
## implementation, which I find unintuitive. Here, if max_miss is set to 1, SNPs with 100%
## missing data would theoretically be allowed, while setting max_miss to 0 will only
## allow SNPs without any missing data
##
## The user can supply a file listing names of samples to keep, one per line. This can
## contain two columns, which is useful if you want to list samples to keep and in the future
## give them new names in column two using the script rename_annotate_split_vcf.sh. Only the
## first column will be used in this script. SAMPLES ARE ALWAYS SORTED in the output,
## regardless of whether sample subsetting is performed.
##
## For depth filtering, this script uses the DP value in the INFO column. This
## value is the total depth of reads that passed the quality control parameters
## defined during SNP calling, (e.g. minimum mapq value). FORMAT DP values, if
## present, are for raw read counts per sample.
##
## snpgap will remove SNPs within n bases of an indel (or overlapping)
## indelgap will thin clusters of indels within n bases of each other, to only retain one
## This script will only retain biallelic variants.
##
## BCFTools is designed to work with .bcf files. If .vcf files are supplied, it
## first converts them to .bcf format internally. Therefore, using a .bcf file
## as input is faster, though BCF files may be larger than .vcf.gz files.
##
## The script will output to either a .bcf file or .vcf.gz depending on the supplied
## file extension
##
## The script will then create summary statistics using bcftools stats, and will
## attempt to create plots of these statistics using the plot-vcfstats script that
## is bundled with bcftools. Python 3 and matplotlib are required for this step.
## If they aren't installed, an error will be printed, but it will not otherwise
## effect the VCF/BCF file output by the script.
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


#### User-defined constants ####

## Note that SNP depth and proportion of missing data are highly correlated

vcf_in="/lustre/project/genolabswheatphg/US_excap/v1_variants/brian_raw_v1_variants/v1_raw_variants.bcf"
out_dir="/lustre/project/genolabswheatphg/US_excap/v1_variants/brian_miss08_maf003_variants/v1_miss08_maf003_variants.bcf"
#regions_bed=""
#vcf_in="/lustre/project/genolabswheatphg/US_excap/v1_variants/1A_filt_test/1A_v1_raw_variants.bcf"
#vcf_out="/lustre/project/genolabswheatphg/US_excap/v1_variants/1A_filt_test/new_code_1A_v1_filt_variants.bcf"
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

array_ind=$1
regions_bed=$2

## Get first letter of true/false het2miss option
## Set some constants depending on het2miss
het2miss="${het2miss:0:1}"
if [[ "$het2miss" == [Tt] ]]; then
	max_het="NA"
	het_string="./."
else
	het_string="1/0"
fi


if [[ $array_ind -eq 0 || $array_ind -eq 1 ]]; then
    mkdir -p "$out_dir"
    mkdir -p "${out_dir}/temp_files"

    ## Echo input parameters to output dir
    echo -e "Input VCF\t${vcf_in}" > "${out_dir}/filt_params.txt"
    echo -e "Regions .bed file\t${regions_bed}" >> "${out_dir}/filt_params.txt"
    echo -e "Taxa subset list\t${samp_file}" >> "${out_dir}/filt_params.txt"
    echo -e "Minimum MAF\t${min_maf}" >> "${out_dir}/filt_params.txt"
    echo -e "Max missing proportion\t${max_miss}" >> "${out_dir}/filt_params.txt"
    echo -e "Max het. proportion\t${max_het}" >> "${out_dir}/filt_params.txt"
    echo -e "Min. average depth\t${min_dp}" >> "${out_dir}/filt_params.txt"
    echo -e "Max average depth\t${max_dp}" >> "${out_dir}/filt_params.txt"
    echo -e "Heterozygous calls converted to missing data?\t${het2miss}" >> "${out_dir}/filt_params.txt"
    echo -e "SNP overlap gap\t${snpgap}" >> "${out_dir}/filt_params.txt"
    echo -e "Indel overlap gap\t${indelgap}" >> "${out_dir}/filt_params.txt"

    ## If samp_file exists, use to subset samples
    ## Otherwise retain all samples present in file
    if [[ -f "$samp_file" ]]; then
        cut -f 1 "$samp_file" | sort > "${out_dir}/temp_files/samp_file.txt"
    else
        bcftools query --list-samples "$vcf_in" | sort > "${out_dir}/temp_files/samp_file.txt"
    fi
else
    sleep 5s
fi

## If the supplied array_ind is 0, we copy full regions .bed file to temp dir.
## Otherwise extract one region from .bed file
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



## Get base name and extension of output file; set the output format
#ext="${vcf_out#*.}"
#base="${vcf_out%%.*}"
#if [[ "$ext" == "bcf" ]]; then
#    out_fmt="b"
#elif [[ "$ext" == "gz" ]]; then
#    ext="vcf.gz"
#    out_fmt="z"
#else
#    echo "ERROR - Please supply either a .bcf or .vcf.gz file for output path"
#    exit 1;
#fi



## Create output directory (if necessary) and temp directory
#out_dir=$(dirname "${vcf_out}")
#mkdir -p "${out_dir}"
#temp_dir="$(mktemp -d -p "${out_dir}")"
#if [[ ! "${out_dir}/temp_files" || ! -d "${out_dir}/temp_files" ]]; then
#    echo "ERROR - Could not create temporary directory"
#    exit 1;
#fi


## Echo input parameters to output dir
#echo -e "Input VCF\t${vcf_in}" > "${base}_filt_params.txt"
#echo -e "Output VCF\t${vcf_out}" >> "${base}_filt_params.txt"
#echo -e "Taxa subset list\t${samp_file}" >> "${base}_filt_params.txt"
#echo -e "Minimum MAF\t${min_maf}" >> "${base}_filt_params.txt"
#echo -e "Max missing proportion\t${max_miss}" >> "${base}_filt_params.txt"
#echo -e "Max het. proportion\t${max_het}" >> "${base}_filt_params.txt"
#echo -e "Min. average depth\t${min_dp}" >> "${base}_filt_params.txt"
#echo -e "Max average depth\t${max_dp}" >> "${base}_filt_params.txt"
#echo -e "Unaligned contigs removed?\t${remove_unal}" >> "${base}_filt_params.txt"
#echo -e "Heterozygous calls converted to missing data?\t${het2miss}" >> "${base}_filt_params.txt"
#echo -e "SNP overlap gap\t${snpgap}" >> "${base}_filt_params.txt"
#echo -e "Indel overlap gap\t${indelgap}" >> "${base}_filt_params.txt"



#up_freq=$(echo "1 - $min_maf" | bc -l)

## Big copy-pasted if-else block for different actions depending on remove_unal and het2miss
## Using het2miss will leave you with some SNPs that only have one allele, so some extra steps
## Are necessary to get rid of these
echo "Filtering VCF..."
echo

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
