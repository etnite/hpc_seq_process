#!/bin/bash
set -e

## Combine single sample gVCF files
##
## This script combines single-sample gVCF files contained in an input directory
## into one multi-sample gVCF using GATK's CombineGVCFs tool.
##
## Individual single-sample gVCF files must have ".g.vcf" within their names.
##
## NOTE: Older versions of GATK cannot handle chromosomes as large as those found
## in wheat. Here, a Conda environment named "gatk" is used, as the GATK version
## installed on Ceres is not recent enough. This environment must contain the
## gatk4 package. To make it on Ceres, type:
##
##   module load miniconda
##   conda create --name gatk gatk4
##
## And go through the prompts. In the future, this will likely be unneccessary.
##
## ALSO NOTE: GATK has a new workflow for combining gVCFs - GenomicsDBImport
## This is generally faster/better for downstream variant calling, but for just creating
## a multiple-sample gVCF it requires two steps.
##
## WARNING: GATK is memory-hungry, and will easily overrun
## individual core memory limits if not reigned in. Memory limits should be
## set in the calls to gatk below - for instance, setting:
##
##   gatk --java-options "-Xmx6g" CombineGVCFs ...
##
## will limit memory usage to 6GB.
################################################################################


#### SLURM job control #### 

#SBATCH --job-name="comb-gvcfs" #name of the job submitted
#SBATCH --partition=short #name of the queue you are submitting job to
#SBATCH --nodes=1 #Number of nodes
  ##SBATCH --ntasks=1  #Number of overall tasks - overrides tasks per node
#SBATCH --ntasks-per-node=1 #number of cores/tasks
#SBATCH --time=00:05:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=jane.doe@isp.com #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### User-defined constants ####

in_dir="/project/genolabswheatphg/gvcfs/SRW_single_samp_bw2_excap_GBS_mq20"
ref_gen="/project/genolabswheatphg/v1_refseq/whole_chroms/Triticum_aestivum.IWGSC.dna.toplevel.fa"
out_gvcf="/project/genolabswheatphg/gvcfs/SRW_bw2_excap_GBS_mq20_multisamp.g.vcf.gz"
ncores=1


#### Executable ####

module load parallel
module load samtools
module load bcftools
module load miniconda
source activate gatk

echo
echo "Start time:"
date

## Create output directory, cd to it, create subdirectory for single-chrom gVCFs
out_dir=$(dirname "${out_gvcf}")
mkdir -p "${out_dir}"
cd "${out_dir}"
mkdir -p "${out_dir}"/chrom_gvcfs

## Create list of input gVCF files
printf '%s\n' "${in_dir}"/*.g.vcf.gz > "${out_dir}"/in_gvcfs.list

## Store list of chroms into array
if [[ ! -f "${ref_gen}".fai ]]; then samtools faidx "${ref_gen}"; fi
cut -f 1 "${ref_gen}".fai > "${out_dir}"/chrom_list.txt

## Run CombineGVCFs on chroms in parallel
cat "${out_dir}"/chrom_list.txt | 
    parallel --dryrun -j $ncores "gatk \
         --java-options '-Xmx3g' CombineGVCFs \
         --reference ${ref_gen} \
         --variants in_gvcfs.list \
         --intervals {} \
         --output chrom_gvcfs/chrom_{}.g.vcf.gz" ::: "${chroms[@]}"

## Create list of resulting single-chrom gVCFs
#printf '%s\n' chrom_gvcfs/*.g.vcf.gz > chrom_gvcfs/chr_gvcfs.list

## Concatenate into single output gVCF
#gatk --java-options "-Xmx3g" GatherVcfs --INPUT chrom_gvcfs/chr_gvcfs.list --OUTPUT "${out_gvcf}"
#bcftools index "${out_gvcf}"

#rm -rf "${chrom_gvcfs}"
source deactivate

echo
echo "End time:"
date
