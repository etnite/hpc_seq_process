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
## snpgap will remove SNPs within n bases of an indel (or overlapping)
## indelgap will thin clusters of indels within n bases of each other, to only retain one
################################################################################


#### SLURM job control #### 

#SBATCH --job-name="filt-vcf" #name of the job submitted
#SBATCH --partition=short #name of the queue you are submitting job to
  ##SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1  #Number of overall tasks - overrides tasks per node
  ##SBATCH --ntasks-per-node=22 #number of cores/tasks
#SBATCH --time=06:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=bpward2@ncsu.edu #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### User-defined constants ####

vcf_in="/project/genolabswheatphg/variants/KS_HRW/raw_vcf/KS_HRW_raw.vcf.gz"
vcf_out="/project/genolabswheatphg/variants/KS_HRW/filt_vcf/KS_HRW_filt.vcf.gz"
taxa_list="none"
min_maf=0.03
max_miss=0.8
max_het=0.1
remove_unal="true"
snpgap=3
indelgap=3


#### Executable #####

module load bcftools
module load vcftools

echo
echo "Start time:"
date

## Grab first letter of remove_unal
remove_unal=${remove_unal:0:1}

out_dir=$(dirname ${vcf_out})
mkdir -p $out_dir
temp_dir="$(mktemp -d -p "$out_dir")"
if [[ ! "$temp_dir" || ! -d "$temp_dir" ]]; then
    echo "Could not create temporary directory"
    exit 1;
fi

## Echo input parameters to output dir
echo "Input VCF: ${vcf_in}" > ${out_dir}/filtering_params.txt
echo "Output VCF: ${vcf_out}" >> ${out_dir}/filtering_params.txt
echo "Taxa subset list: ${taxa_list}" >> ${out_dir}/filtering_params.txt
echo "Minimum MAF: ${min_maf}" >> ${out_dir}/filtering_params.txt
echo "Max missing proportion: ${max_miss}" >> ${out_dir}/filtering_params.txt
echo "Max het. proportion: ${max_het}" >> ${out_dir}/filtering_params.txt
echo "Unaligned contigs removed? ${remove_unal}" >> ${out_dir}/filtering_params.txt
echo "SNP overlap gap: ${snpgap}" >> ${out_dir}/filtering_params.txt
echo "Indel overlap gap: ${indelgap}" >> ${out_dir}/filtering_params.txt

## If taxa_list exists, use to subset samples
## Otherwise retain all samples present in VCF file
if [[ -f $taxa_list ]]; then
	cp $taxa_list $temp_dir/taxa_list.txt
else
	bcftools query -l $vcf_in > $temp_dir/taxa_list.txt
fi

## Invert the max missing parameter
snpmissinv="$(echo "1 - ${max_miss}" | bc)"

## First round of filtering for taxa list, missing data and MAF
echo "Filtering by missing data and MAF..."
echo
if [[ $remove_unal == [Tt] ]]; then
	vcftools --gzvcf $vcf_in \
	         --keep $temp_dir/taxa_list.txt \
		     --max-missing $snpmissinv \
		     --maf $min_maf \
		     --not-chr UN \
		     --recode \
		     --stdout > ${temp_dir}/phase1_filt.vcf
elif [[ $remove_unal == [Ff] ]]; then
	vcftools --gzvcf $vcf_in \
	         --keep $temp_dir/taxa_list.txt \
		     --max-missing $snpmissinv \
		     --maf $min_maf \
		     --recode \
		     --stdout > ${temp_dir}/phase1_filt.vcf
else
	echo "Please supply 'true' or 'false' for remove_unal"
	exit 1;
fi

## Find SNPs that pass the specified het level
## Calculate number of non-missing calls per SNP
echo "Calculating non-missing calls per site..."
echo
bcftools query -f '[\S%CHROM\_%POS\n]' -i 'GT!~".\."' ${temp_dir}/phase1_filt.vcf | sort -T $temp_dir | uniq -c > ${temp_dir}/non_miss.txt

## Calculate number of het calls per SNP
echo "Calculating heterozygosity per site..."
echo
bcftools query -f '[\S%CHROM\_%POS\n]' -i 'GT="het"' ${temp_dir}/phase1_filt.vcf | sort -T $temp_dir | uniq -c > ${temp_dir}/het.txt

## Join the two files together, inserting 0 for missing vals (should only be
## values from the het file that have missing)
## Format is: SNP_name het_count non-missing_count
join -a 2 -e 0 -1 2 -2 2 -o 2.2,1.1,2.1 ${temp_dir}/het.txt ${temp_dir}/non_miss.txt > ${temp_dir}/joined.txt

## Find SNPs passing the het. threshold
awk -v het="${max_het}" '$2/$3 <= het {print $1}' ${temp_dir}/joined.txt > ${temp_dir}/snps_keep.txt

## Second round of filtering for heterozygosity
echo "Filtering by site heterozygosity..."
echo
vcftools --vcf ${temp_dir}/phase1_filt.vcf \
         --snps ${temp_dir}/snps_keep.txt \
         --recode \
         --stdout | 
         bcftools filter -g $snpgap -G $indelgap -Oz - -o $vcf_out
bcftools index -c $vcf_out


## Generate summary stats using TASSEL
#echo "Generating summary statistics with TASSEL"
#$TASSEL_PL -vcf "${vcf_out}" \
#           -GenotypeSummaryPlugin \
#           -endPlugin \
#           -export "${out_dir}"/summary

rm -rf $temp_dir
#source deactivate

echo
echo "End time:"
date
