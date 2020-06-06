#!/bin/bash
set -e


## Concatenate and sort VCF files
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script concatenates and sorts a set of VCF or BCF files, producing a
## single compressed BCF file as output.
##
## In addition, it adds SNP IDs formatted S<chromosome>_<position>, left aligns
## indels, and normalizes reference alleles.
################################################################################


#### SLURM job control #### 

#SBATCH --job-name="cat-sort-vcfs" #name of the job submitted
#SBATCH --partition=mem #name of the queue you are submitting job to
#SBATCH --nodes=1 #Number of nodes
  ##SBATCH --ntasks=28  #Number of overall tasks - overrides tasks per node
#SBATCH --ntasks-per-node=22 #number of cores/tasks
#SBATCH --time=48:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=jane.doe@isp.com #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### User-defined constants ####

vcfs_list="/project/genolabswheatphg/alignments/ERSGGL_SRW_bw2_bams/SRW_merged_excap_GBS_wholechrom_bw2_bams"
out_bcf="/project/genolabswheatphg/variants/SRW/ERSGGL_SRW_merged_excap_GBS_wholechr_bw2.vcf.gz"

## Maximum memory to use, which can be set to, e.g. "3G" for 3 gigabytes
max_memory="10G"


#### Executable ####

module load bcftools

echo
echo "Start time:"
date

out_dir=$(dirname "$out_bcf")
mkdir -p "$out_dir"

## Normalize indels; discard all but one overlapping SNP, all but one overlapping indel
bcftools concat --no-version \
                --file-list "$vcfs_list" \
                --output-type u |
    bcftools sort --temp-dir "$out_dir" \
                  --max-mem "$max_memory" \
                  --output-type u |
    bcftools annotate --set-id +'S%CHROM\_%POS' \
                      --output-type u |
    bcftools norm --fasta-ref "$ref_gen" \
                  --threads $SLURM_NTASKS \
                  --output-type b > "$out_bcf"
 
bcftools index -c "$out_bcf"

echo
echo "End time:"
date
