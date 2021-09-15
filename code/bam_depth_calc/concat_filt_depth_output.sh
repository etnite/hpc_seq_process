#!/bin/bash

## Merge samtools_depth_parallel.sh output
##
## This is just a simple script to turn the output of the associated script
## samtools_depth_parallel.sh into a single, more compact .bed file. The bedtools
## merge tools supports adding in various notations in an additional fourth
## column. This gets a little tricky to interpret though. In the example below,
## we filter our samtools depth output by the median read depth at each position.
## Then we merge together everything separated by less than the length of a read
## (in this case 150bp). In the output .bed file we also calculate the mean of
## column 4, which is the number of reads at that position. Note that this IS NOT
## the mean read depth at a particular position - it is the mean of the read
## depths across the entire interval listed in the output .bed file.
################################################################################


#### SLURM job control parameters ####

## Options with two comment chars are deactivated

#SBATCH --account=guedira_seq_map
#SBATCH --job-name="bed-cat"  #name of the job submitted
#SBATCH --partition=atlas  #name of the queue you are submitting job to
#SBATCH --ntasks=1  #Number of overall tasks - overrides tasks per node
#SBATCH --time=00:05:00 #time allocated for this job hours:mins:seconds
    ##SBATCH --mail-user=ward.1660@osu.edu #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### User-Defined Constants ####

in_dir="/project/guedira_seq_map/Allegro_test/bam_R_depth_calc_mq20_min_mdn4/"
out_bed="/project/guedira_seq_map/Allegro_test/groupA_R_strand_depths_mq20_min_mdn4.bed"
min_mdn=4
merge_dist=150


#### Executable ####

echo
echo "Start time:"
date

shopt -s extglob

module load singularity/3.7.1
module load bedtools

cat "$in_dir"/!(*header).txt |
    awk -v min_med=$min_mdn 'BEGIN{OFS = "\t"} $7 >= min_med {print $1, $2, $3, $4}' |
    bedtools merge -i stdin -d $merge_dist -c 4 -o mean > "$out_bed"

echo
echo "End time:"
date
