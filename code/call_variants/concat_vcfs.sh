#!/bin/bash


## Concatenate VCF files
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script is identical to concat_sort_vcfs.sh except that it does not perform
## any sorting step - this is intended for situations where the names of the VCF
## are sufficient for proper ordering during the concatenation step.
##
## In addition, it adds SNP IDs formatted S<chromosome>_<position>, left aligns
## indels, and normalizes reference alleles.
##
## NOTE: reference genome must be indexed with samtools faidx
################################################################################


#### SLURM job control #### 

#SBATCH -A guedira_seq_map
#SBATCH --job-name="cat-vcfs" #name of the job submitted
#SBATCH --partition=atlas #name of the queue you are submitting job to
#SBATCH --nodes=1 #Number of nodes
  ##SBATCH --ntasks=28  #Number of overall tasks - overrides tasks per node
#SBATCH --ntasks-per-node=10 #number of cores/tasks
#SBATCH --time=01:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=jane.doe@isp.com #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### User-defined constants ####

## Path to text file listing input VCFs or BCFs, one per line
## Path to reference genome .fasta file
## Path to write output BCF file to
vcfs_list="/home/brian.ward/samp_and_file_lists/groupB_Rstrand_region_bcfs_list.txt"
ref_gen="/project/guedira_seq_map/ref_genomes/v1_refseq_w_KIMs/CSv1_refseq_w_KIMs.fa"
out_bcf="/project/guedira_seq_map/Allegro_test/groupB_mq20_raw_VCFs/R_strand/groupB_mq20_R_strand_raw.bcf"


#### Executable ####

## Using conda here to utilize bcftools plugins and plot-vcfstats
module load singularity/3.7.1
module load bcftools
#module load miniconda
#source activate htslib

echo
echo "Start time:"
date

out_dir=$(dirname "$out_bcf")
mkdir -p "$out_dir"
base="${out_bcf%%.*}"

## Normalize indels; discard all but one overlapping SNP, all but one overlapping indel
bcftools concat --no-version \
    --file-list "$vcfs_list" \
    --output-type u |
bcftools annotate --no-version \
    --set-id +'S%CHROM\_%POS' \
    --output-type u |
bcftools norm --no-version \
    --fasta-ref "$ref_gen" \
    --threads $SLURM_NTASKS \
    --output-type b > "$out_bcf"
 
bcftools index -c "$out_bcf"

## Generate summary stats
bcftools stats "${out_bcf}" > "${base}_stats.txt"

## Create plots from summary stats
mkdir "${base}_plots"
plot-vcfstats --prefix "${base}_plots" --no-PDF "${base}_stats.txt"

echo
echo "End time:"
date
