#!/bin/bash
set -e

## VCF Beagle Imputation Without Reference Panel
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script takes an input BCF or VCF file, runs imputation on it using Beagle
## and outputs a bgzipped VCF file. Note that the output is always a .vcf.gz
## regardless of the input format.
##
## This script expects the path to the Beagle .jar archive to be aliased as
## "beajar". This can be set either in the user's ~/.bashrc or ~/.bash_profile
## file.
##
## Currently this script also expects a genetic map (I use the Synthetic x Opata
## GBS map of Gutierrez-Gonzalez et al. [DOI 10.1007/s00122-004-1740-7]). The script
## can accept a Plink-formatted genetic map file for the map_file parameter.
## However, this can be set to any string that isn't a valid file to run Beagle
## without a genetic map.
##
## The script is configured to run entire wheat chromosomes at once. To disable
## this, either modify the window argument in the lines starting with java -jar, or
## else comment them out.
##
## The output is a bgzipped VCF file in the same directory as the input file, with
## "_imp" inserted before the file extension. 
################################################################################


#### SLURM job control #### 

#SBATCH --job-name="beagle-imp" #name of the job submitted
#SBATCH --partition=short #name of the queue you are submitting job to
#SBATCH --nodes=1 #Number of nodes
    ##SBATCH --ntasks=1  #Number of overall tasks - overrides tasks per node
#SBATCH --ntasks-per-node=20 #number of cores/tasks
#SBATCH --time=06:00:00 time allocated for this job hours:mins:seconds
#SBATCH --mail-user=bpward2@ncsu.edu #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### Read input parameters ####

vcf_in=$1
map_file=$2
n_threads=$3


#### Executable ####

#module load java
#module load bcftools

echo
echo "Start time:"
date

## Get base name and extension of input file; set output file name
ext="${vcf_in#*.}"
base="${vcf_in%%.*}"
base="${vcf_in%%.*}"
vcf_out="${base}_imp"

## If working with a BCF file, we first have to convert to VCF
if [[ "$ext" == "bcf" ]]; then
    bcftools view "$vcf_in" -Oz -o "${base}.vcf.gz"
    vcf_in="${base}.vcf.gz"
elif [[ "$ext" == "gz" || "$ext" == "vcf.gz" ]]; then
    :
else
    echo "ERROR - Please supply either a .bcf or .vcf.gz file for input path"
    exit 1;
fi

## Run Beagle with or without map file
echo
if [[ -f "$map_file" ]]; then
    echo "Imputing with map file"    
    java -jar $beajar gt="$vcf_in" out="$vcf_out" map="$map_file" nthreads=$n_threads window=205
else
    echo "Imputing without map file"
    java -jar $beajar gt="$vcf_in" out="$vcf_out" nthreads=$n_threads window=205
fi

## Update the header in the output VCF to include contig info
bcftools view --header-only "${vcf_out}.vcf.gz" | grep "^##" > "${base}_new_header.txt"
bcftools view --header-only "$vcf_in" | grep "^##contig" >> "${base}_new_header.txt"
bcftools view --header-only "${vcf_out}".vcf.gz | grep "^#CHROM" >> "${base}_new_header.txt"
bcftools reheader "${vcf_out}.vcf.gz" --header "${base}_new_header.txt" --output "${base}_temp.vcf.gz"
mv "${base}_temp.vcf.gz" "${vcf_out}.vcf.gz"

rm "${base}_new_header.txt"

bcftools index "${vcf_out}.vcf.gz"

echo
echo "End time:"
date
