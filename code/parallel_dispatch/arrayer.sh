#!/bin/bash

## Parallelize script with SLURM job arrays
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script is a "universal parallelizer" - its intended function is to
## enable the parallel execution of another script WHEN DIFFERENT PROCESSES
## CAN BE RUN INDEPENDENTLY, such as when running the same process on different
## samples, as is common in bioinformatics.
##
## This script uses an array of integers (set in the line starting with
## "#SBATCH --array=") to perform parallel dispatch for the script it is calling.
## These integers are typically used to subset out individual lines of a file
## that is iterated over. For instance, setting "#SBATCH --array=1-10" will
## subset out lines 1 - 10 of the file being iterated over, and supply them to
## ten different instances of the script being called. In this case, each
## instance of the script that is run would receive one of the first ten lines
## of the iterated file as a positional input argument. Typically, the file that
## is iterated over lists samples, filenames, or genomic regions.
##
## The user can optionally specify the file to iterate over as iter_file. If this
## is set to any string that is not a valid filename, then it assumes that the
## path to the file to iterate over is supplied within the script being called.
################################################################################


#### SLURM job control parameters ####

## Options with two comment chars are deactivated

#SBATCH --job-name="arrayer"  #name of the job submitted
#SBATCH --partition=short  #name of the queue you are submitting job to
#SBATCH --array=1-59%5  #array range - can choose number simultaneous jobs with %, e.g. --array=1-12%4
#SBATCH --nodes=1 #Number of nodes
  ##SBATCH --ntasks=28  #Number of overall tasks - overrides tasks per node
#SBATCH --ntasks-per-node=10 #number of cores/tasks
#SBATCH --time=06:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=jane.doe@isp.com #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### User-Defined Constants ####

## To specify file to iterate through here, supply its path for iter_file.
## Otherwise set iter_file to a string that is not a valid file name, like "none"
## or "nothing". In this case the file to iterate over should be set in the
## script that is being called.
script="alignment/bowtie2_align_parallel.sh"
iter_file="nothing"


#### Executable ####

echo
echo "${script}"
echo "Start time:"
date

if [[ -f "$iter_file" ]]; then
    bash "$script" "$iter_file" $SLURM_ARRAY_TASK_ID
else
    bash "$script" $SLURM_ARRAY_TASK_ID
fi

echo
echo "End time:"
date
