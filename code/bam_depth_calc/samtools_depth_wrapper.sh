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

#SBATCH --job-name="depth-stats" #name of the job submitted
#SBATCH --partition=short #name of the queue you are submitting job to
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1  #Number of overall tasks - overrides tasks per node
#SBATCH --ntasks-per-node=22 #number of cores/tasks
#SBATCH --time=24:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=jane.doe@isp.com #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### User-Defined Constants ####

bam_list_file="gbs_test_bamlist.txt"
max_dep=200
min_mean=0
min_median=5
max_dist=90
#min_size=50
#out_file="depth_test_out.txt"


#### Executable ####

#module load samtools
#module load miniconda
#source activate py38

samtools depth -m $max_dep -f "$bam_list_file" |
    ./samtools_depth_parse.jl -u $min_mean -m $min_median -d $max_dist |
    head -n 1000 > test_out.bed


#    ./samtools_depth_summstats.py |
#    gzip -c > "$out_file"

#source deactivate
