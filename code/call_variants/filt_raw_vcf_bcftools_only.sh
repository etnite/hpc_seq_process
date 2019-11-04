#!/bin/bash
set -e
#source /home/gbg_lab_admin/miniconda3/bin/activate bwa_align_call


## Filter VCF file - single threaded
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
## as input is faster, though BCF files are often larger than .vcf.gz files.
################################################################################


#### SLURM job control #### 

#SBATCH --job-name="filt-vcf" #name of the job submitted
#SBATCH --partition=short #name of the queue you are submitting job to
  ##SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1  #Number of overall tasks - overrides tasks per node
  ##SBATCH --ntasks-per-node=22 #number of cores/tasks
#SBATCH --time=06:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=jane.doe@isp.com #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### User-defined constants ####

## Note that SNP depth and proportion of missing data are highly correlated

vcf_in="/project/genolabswheatphg/variants/SRW/ERSGGL_SRW_merged_excap_GBS_wholechr_bw2.vcf.gz"
vcf_out="/project/genolabswheatphg/variants/SRW/ERSGGL_SRW_merged_excap_GBS_wholechr_bw2_filt08.vcf.gz"
taxa_list="none"
min_maf=0.03
max_miss=0.8
max_het=0.1
min_dp=0
max_dp=1e6
remove_unal="false"
snpgap=3
indelgap=3


#### Executable #####

module load bcftools
#module load vcftools

echo
echo "Start time:"
date

## Grab first letter of remove_unal
remove_unal=${remove_unal:0:1}

out_dir=$(dirname "${vcf_out}")
mkdir -p "${out_dir}"
temp_dir="$(mktemp -d -p "${out_dir}")"
if [[ ! "${temp_dir}" || ! -d "${temp_dir}" ]]; then
    echo "Could not create temporary directory"
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
echo -e "SNP overlap gap\t${snpgap}" >> "${out_dir}"/filtering_params.txt
echo -e "Indel overlap gap\t${indelgap}" >> "${out_dir}"/filtering_params.txt

## If taxa_list exists, use to subset samples
## Otherwise retain all samples present in VCF file
if [[ -f $taxa_list ]]; then
	cp $taxa_list "${temp_dir}"/taxa_list.txt
else
	bcftools query --list-samples $vcf_in > "${temp_dir}"/taxa_list.txt
fi

## Some real bcftools power-user stuff here
echo "Filtering VCF..."
echo
if [[ $remove_unal == [Tt] ]]; then
    bcftools view "${vcf_in}" \
        --samples-file "${temp_dir}"/taxa_list.txt \
        --output-type u |
    bcftools view - \
        --targets ^UN \
        --exclude "F_MISSING > ${max_miss} || MAF < ${min_maf} || INFO/DP < ${min_dp} || INFO/DP > ${max_dp} || (COUNT(GT=\"het\") / COUNT(GT!~\"\.\")) > ${max_het} " \
        --output-type u |
    bcftools filter --SnpGap $snpgap \
        --IndelGap $indelgap \
        --output-type z \
        --output "${vcf_out}"
elif [[ $remove_unal == [Ff] ]]; then
	bcftools view "${vcf_in}" \
        --samples-file "${temp_dir}"/taxa_list.txt \
        --output-type u |
    bcftools view - \
        --exclude "F_MISSING > ${max_miss} || MAF < ${min_maf} || INFO/DP < ${min_dp} || INFO/DP > ${max_dp} || (COUNT(GT=\"het\") / COUNT(GT!~\"\.\")) > ${max_het} " \
        --output-type u |
    bcftools filter --SnpGap $snpgap \
        --IndelGap $indelgap \
        --output-type z \
        --output "${vcf_out}"
else
	echo "Please supply 'true' or 'false' for remove_unal"
	exit 1;
fi

bcftools index -c "${vcf_out}"


## Generate summary stats using TASSEL
#echo "Generating summary statistics with TASSEL..."
#$TASSEL_PL -vcf "${vcf_out}" \
#           -GenotypeSummaryPlugin \
#           -endPlugin \
#           -export "${out_dir}"/summary

rm -rf $temp_dir
#source deactivate

echo
echo "End time:"
date
