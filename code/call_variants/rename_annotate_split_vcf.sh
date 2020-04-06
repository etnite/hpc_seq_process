#!/bin/bash
set -e
#source /home/gbg_lab_admin/miniconda3/bin/activate bwa_align_call


## Rename samples, annotate variants, split SNPs and indels
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## INPUT FILES:
##
##   1. VCF file
##   2. reference .fasta file (must be indexed with samtools)
##   3. .gff3 annotations file
##   4. (optional) two-column tab-delimited file, with current sample names in
##      first column, and corresponding desired sample names in second column. 
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
#SBATCH --time=02:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=jane.doe@isp.com #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### User-defined constants ####

vcf_in="/project/genolabswheatphg/variants/KS_HRW/filt_vcf/KS_HRW_filt.vcf.gz"
ref="/project/genolabswheatphg/v1_refseq/whole_chroms/Triticum_aestivum.IWGSC.dna.toplevel.fa"
gff="/project/genolabswheatphg/v1_refseq/whole_chroms/Triticum_aestivum.IWGSC.44.gff3.gz"
lookup_file="none"


#### Executable ####

module load bcftools

echo
echo "Start time:"
date

out_dir=$(dirname ${vcf_in})
base=$(echo ${vcf_in} | sed 's/.vcf.gz$//')

## If sample name lookup file is supplied, then replace sample names
## and call consequences. Otherwise, skip sample renaming step
if [[ -f $lookup_file ]]; then
	bcftools reheader --samples $lookup_file $vcf_in |
        bcftools csq -f $ref -g $gff -p a -Oz -o ${base}_csq.vcf.gz
else
    bcftools csq -f $ref -g $gff -p a -Oz $vcf_in -o ${base}_csq.vcf.gz
fi

bcftools index -c ${base}_csq.vcf.gz

## Create SNP-only and indel-only files
base=$(echo ${vcf_in} | sed 's/.vcf.gz$//')
bcftools view --type snps -Oz -o ${base}_csq_snps_only.vcf.gz ${base}_csq.vcf.gz
bcftools index -c ${base}_csq_snps_only.vcf.gz
bcftools view --type indels -Oz -o ${base}_csq_indels_only.vcf.gz ${base}_csq.vcf.gz
bcftools index -c ${base}_csq_indels_only.vcf.gz

#source deactivate

echo
echo "End time:"
date
