#!/bin/bash

## Wrapper script for samtools_depth_summstats.py
##
## This script handles input and output streaming for
## samtools_depth_summstats.py
##
## samtools_depth_summstats.py requires numpy, so here miniconda is
## used to load the environment "py38" which contains the dependency
####################################################################


#### SLURM job control #### 

#SBATCH --job-name="jul-depth-stats" #name of the job submitted
#SBATCH --partition=short #name of the queue you are submitting job to
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1  #Number of overall tasks - overrides tasks per node
    ##SBATCH --ntasks-per-node=22 #number of cores/tasks
#SBATCH --time=48:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=bpward2@ncsu.edu #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### User-Defined Constants ####

## Can provide a file listing bam files, one per line
## Or supply a string that isn't a filename for bam_list_file
## And supply a directory to look for .bam files
bam_list_file="none"
bam_dir="/project/genolabswheatphg/alignments/ERSGGL_SRW_bw2_bams/SRW_merged_excap_GBS_wholechrom_bw2_bams_mq20_filt"

max_dep=200
min_qual=20
min_mean=4.31
min_median=0
max_dist=150
#min_size=50

out_file="/project/genolabswheatphg/SRW_depth_test/julia_test/julia_1A_depth_regions.bed"


#### Executable ####

module load samtools

echo
echo "Start time:"
date

## Create list of .bam files
out_dir=$(dirname "$out_file")
mkdir -p "$out_dir"
if [[ -f "$bam_list_file" ]]; then
    cp "$bam_list_file" "$out_dir"/temp_bam_list.txt
else
    realpath "$bam_dir"/*.bam > "$out_dir"/temp_bam_list.txt
fi

## Run depth calculation
samtools depth -r 1A -m $max_dep -Q $min_qual -f "$out_dir"/temp_bam_list.txt |
    ./samtools_depth_parse.jl -u $min_mean -m $min_median -d $max_dist > "$out_file"


#    ./samtools_depth_summstats.py |
#    gzip -c > "$out_file"

rm "$out_dir"/temp_bam_list.txt

echo
echo "End time:"
date
