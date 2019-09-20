#!/bin/bash

## Parallelize script with SLURM job arrays
##
## This script is a "universal parallelizer" - it's intended function is to
## enable the parallel execution of another script WHEN DIFFERENT PROCESSES
## CAN BE RUN INDEPENDENTLY, such as when running the same process on different
## samples, as is common in bioinformatics.
##
## This script must read in an array to use as the parallel iterator. This is
## typically a text file containing a list of some kind with one entry per line.
## Each entry in this list is then fed into whatever script is being called as a
## positional argument. SLURM then handles the parallel execution. Additional 
## parameters often must be set in the script that is being called.
################################################################################


#### SLURM job control parameters ####

## Options with two comment chars are deactivated

#SBATCH --job-name="bw2-align"  #name of the job submitted
#SBATCH --partition=short  #name of the queue you are submitting job to
#SBATCH --array=1  #array range - can choose number simultaneous jobs with %, e.g. --array=1-12%4
#SBATCH --nodes=1 #Number of nodes
  ##SBATCH --ntasks=28  #Number of overall tasks - overrides tasks per node
#SBATCH --ntasks-per-node=20 #number of cores/tasks
#SBATCH --time=06:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=bpward2@ncsu.edu #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### User-Defined Constants ####

script="/home/brian.ward/repos/wheat_phg/code/alignment/bowtie2_SE_align_parallel.sh"


#### Executable ####

echo
echo "${script}"
echo "Start time:"
date

bash $script $SLURM_ARRAY_TASK_ID

echo
echo "End time:"
date
