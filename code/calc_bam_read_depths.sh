#!/bin/bash
set -e

## Calculate read depth for BAM files
##
## This generates a file listing read depth at all positions within a genome,
## using all BAM files within a user-defined directory. The bam files are
## expected to be in the format input-directory/samplename.bam.
##
## The output is a single gzipped, tab-delimited file, with chromosome as the
## first column, position as second, and remaining columns read depths across
## all included samples. Positions with read depths of 0 across all samples are
## excluded
################################################################################


#### SLURM job control #### 

#SBATCH --job-name="read-depths" #name of the job submitted
#SBATCH --partition=short #name of the queue you are submitting job to
#SBATCH --ntasks=28  #Number of overall tasks - overrides tasks per node
#SBATCH --time=8:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=bpward2@ncsu.edu #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### User-defined constants ####

in_dir="my_bam_files"
out_file="depth_test.tsv.gz"


#### Executable ####

module load samtools

echo
echo "Start time:"
date

## Create output directory
out_dir=$(dirname "${out_vcf}")
mkdir -p "${out_dir}"

## Sanity check on output filename
ext="${out_file##*.}"
if [[ $ext != "gz" ]]; then
    echo -e "\nOutput filename must have .gz extension"
    exit 1;
fi

## Create list of .bam files
printf '%s\n' "${in_dir}"/*.bam > "${out_dir}"/bam_list.txt

## Create sample list and header line for output file
echo "chr" > "${out_dir}"/samp_list.txt
echo "pos" >> "${out_dir}"/samp_list.txt
sed 's/.*\///' bam_list.txt | sed 's/.bam/_DP/' >> "${out_dir}"/samp_list.txt
cat "${out_dir}"/samp_list.txt | tr '\n' '\t' | gzip -c > "${out_file}"

## Run the depth calculation
samtools depth -f "${out_dir}"/bam_list.txt | gzip -c >> "${out_file}"

echo
echo "End time:"
date
