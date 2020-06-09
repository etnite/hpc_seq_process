#!/bin/bash
set -e


## Filter VCF file - single threaded
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## NOTE: The max_miss parameter works in an opposite manner of the VCFTools
## implementation, which I find unintuitive. Here, if max_miss is set to 1, SNPs with 100%
## missing data would theoretically be allowed, while setting max_miss to 0 will only
## allow SNPs without any missing data
##
## For depth filtering, this script uses the DP value in the INFO column. This
## value is the total depth of reads that passed the quality control parameters
## defined during SNP calling, (e.g. minimum mapq value). FORMAT DP values, if
## present, are for raw read counts per sample.
##
## snpgap will remove SNPs within n bases of an indel (or overlapping)
## indelgap will thin clusters of indels within n bases of each other, to only retain one
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
## effect the output of the script.
################################################################################


#### SLURM job control #### 

#SBATCH --job-name="filt-vcf" #name of the job submitted
#SBATCH --partition=short #name of the queue you are submitting job to
  ##SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1  #Number of overall tasks - overrides tasks per node
  ##SBATCH --ntasks-per-node=22 #number of cores/tasks
#SBATCH --time=5:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=bpward2@ncsu.edu #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### User-defined constants ####

## Note that SNP depth and proportion of missing data are highly correlated

vcf_in="/project/guedira_seq_map/brian/US_excap/v1_variants/US_excap_raw_variants.bcf"
vcf_out="/project/guedira_seq_map/brian/US_excap/v1_variants/US_excap_filt08_variants.bcf"
taxa_list="none"
min_maf=0.05
max_miss=0.8
max_het=1
min_dp=0
max_dp=1e6
remove_unal="false"
het2miss="true"
snpgap=3
indelgap=3


#### Executable #####

## Using conda here to utilize bcftools plugins and plot-vcfstats
#module load bcftools
module load miniconda
source activate htslib

echo
echo "Start time:"
date

## Get base name and extension of output file; set the output format
ext="${vcf_out#*.}"
base="${vcf_out%%.*}"
if [[ "$ext" == "bcf" ]]; then
    out_fmt="b"
elif [[ "$ext" == "vcf.gz" ]]; then
    out_fmt="z"
else
    echo "ERROR - Please supply either a .bcf or .vcf.gz file for output path"
    exit 1;
fi

## Get first letter of true/false options
remove_unal="${remove_unal:0:1}"
het2miss="${het2miss:0:1}"
if [[ "$het2miss" == [Tt] ]]; then max_het="NA"; fi

## Create output directory (if necessary) and temp directory
out_dir=$(dirname "${vcf_out}")
mkdir -p "${out_dir}"
temp_dir="$(mktemp -d -p "${out_dir}")"
if [[ ! "${temp_dir}" || ! -d "${temp_dir}" ]]; then
    echo "ERROR - Could not create temporary directory"
    exit 1;
fi

## Echo input parameters to output dir
echo -e "Input VCF\t${vcf_in}" > "${out_dir}"/filtering_params.txt
echo -e "Output VCF\t${vcf_out}" >> "${out_dir}"/filtering_params.txt
echo -e "Taxa subset list\t${taxa_list}" >> "${out_dir}"/filtering_params.txt
echo -e "Minimum MAF\t${min_maf}" >> "${out_dir}"/filtering_params.txt
echo -e "Max missing proportion\t${max_miss}" >> "${out_dir}"/filtering_params.txt
echo -e "Max het. proportion\t${max_het}" >> "${out_dir}"/filtering_params.txt
echo -e "Min. average depth\t${min_dp}" >> "${out_dir}"/filtering_params.txt
echo -e "Max average depth\t${max_dp}" >> "${out_dir}"/filtering_params.txt
echo -e "Unaligned contigs removed?\t${remove_unal}" >> "${out_dir}"/filtering_params.txt
echo -e "Heterozygous calls converted to missing data?\t${het2miss}" >> "${out_dir}"/filtering_params.txt
echo -e "SNP overlap gap\t${snpgap}" >> "${out_dir}"/filtering_params.txt
echo -e "Indel overlap gap\t${indelgap}" >> "${out_dir}"/filtering_params.txt

## If taxa_list exists, use to subset samples
## Otherwise retain all samples present in VCF file
if [[ -f $taxa_list ]]; then
    cp $taxa_list "${temp_dir}"/taxa_list.txt
else
    bcftools query --list-samples $vcf_in > "${temp_dir}"/taxa_list.txt
fi


## Big if-else block for different actions depending on remove_unal and het2miss
echo "Filtering VCF..."
echo
if [[ "$remove_unal" == [Tt] && "$het2miss" == [Tt] ]]; then
    bcftools view "${vcf_in}" \
        --samples-file "${temp_dir}"/taxa_list.txt \
        --output-type u |
    bcftools +setGT --output-type u - -- --target-gt q \
        --include 'GT="het"' \
        --new-gt "./." |
    bcftools view - \
        --targets ^UN,Un \
        --exclude "F_MISSING > ${max_miss} || MAF < ${min_maf} || INFO/DP < ${min_dp} || INFO/DP > ${max_dp}" \
        --output-type u |
    bcftools filter --SnpGap $snpgap \
        --IndelGap $indelgap \
        --output-type "$out_fmt" \
        --output "${vcf_out}"
elif [[ "$remove_unal" == [Tt] && "$het2miss" == [Ff] ]]; then
    bcftools view "${vcf_in}" \
        --samples-file "${temp_dir}"/taxa_list.txt \
        --output-type u |
    bcftools view - \
        --targets ^UN,Un \
        --exclude "F_MISSING > ${max_miss} || MAF < ${min_maf} || INFO/DP < ${min_dp} || INFO/DP > ${max_dp} || (COUNT(GT=\"het\") / COUNT(GT!~\"\.\")) > ${max_het} " \
        --output-type u |
    bcftools filter --SnpGap $snpgap \
        --IndelGap $indelgap \
        --output-type "$out_fmt" \
        --output "${vcf_out}"
elif [[ $remove_unal == [Ff] && "$het2miss" == [Tt] ]]; then
    bcftools view "${vcf_in}" \
        --samples-file "${temp_dir}"/taxa_list.txt \
        --output-type u |
    bcftools +setGT --output-type u - -- --target-gt q \
        --include 'GT="het"' \
        --new-gt "./." |
    bcftools view - \
        --exclude "F_MISSING > ${max_miss} || MAF < ${min_maf} || INFO/DP < ${min_dp} || INFO/DP > ${max_dp}" \
        --output-type u |
    bcftools filter --SnpGap $snpgap \
        --IndelGap $indelgap \
        --output-type "$out_fmt" \
        --output "${vcf_out}"
elif [[ $remove_unal == [Ff] && "$het2miss" == [Ff] ]]; then
	bcftools view "${vcf_in}" \
        --samples-file "${temp_dir}"/taxa_list.txt \
        --output-type u |
    bcftools view - \
        --exclude "F_MISSING > ${max_miss} || MAF < ${min_maf} || INFO/DP < ${min_dp} || INFO/DP > ${max_dp} || (COUNT(GT=\"het\") / COUNT(GT!~\"\.\")) > ${max_het} " \
        --output-type u |
    bcftools filter --SnpGap $snpgap \
        --IndelGap $indelgap \
        --output-type "$out_fmt" \
        --output "${vcf_out}"
else
	echo "ERROR - Please supply 'true' or 'false' for remove_unal and het2miss"
	exit 1;
fi

bcftools index -c "${vcf_out}"

## Generate summary stats
bcftools stats "${vcf_out}" > "${base}_stats.txt"

## Create plots from summary stats
mkdir "${base}_plots"
plot-vcfstats --prefix "${base}_plots" --no-PDF "${base}_stats.txt"

## Generate summary stats using TASSEL
#echo "Generating summary statistics with TASSEL..."
#$TASSEL_PL -vcf "${vcf_out}" \
#           -GenotypeSummaryPlugin \
#           -endPlugin \
#           -export "${out_dir}"/summary

rm -rf $temp_dir


echo
echo "End time:"
date
