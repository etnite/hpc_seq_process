#!/bin/bash

## Parallelize script with GNU Parallel
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script is a "universal parallelizer" - it's intended function is to
## enable the parallel execution of another script WHEN DIFFERENT PROCESSES
## CAN BE RUN INDEPENDENTLY, such as when running the same process on different
## samples, as is common in bioinformatics.
##
## The script arrayer.sh offers finer control on HPC clusters - e.g. the ability
## to define the number of cores used per process, and to control the number of
## simultaneously running jobs.
##
## This script defines an array of integers, from one to a user-defined maximum.
## Typically, each integer will then be used to subset a single line from a file
## that might list samples, files, genomic regions, etc.
##
## NOTES control parameters for jobs:
##
## Jobs that are dispatched with GNU parallel on a HPC cluster WILL NOT inherit
## the parent script's SLURM constants (unlike using SLURM arrays with arrayer.sh).
## Therefore, constants such as number of tasks per job or time limit must be set within the
## child script. Note that all jobs will be submitted simultaneously (unlike
## arrayer.sh, which allows for setting the number of simultaneously running jobs).
##
## NOTES on using on single machine vs. cluster:
##
## If running on a single machine, the line "module load parallel" should be
## commented out (GNU parallel must be installed on the machine). All SLURM job
## control lines will be treated as comments and ignored.
##
## The user can optionally specify the file to iterate over as iter_file. If this
## is set to any string that is not a valid filename, then it assumes that the
## path to the file to iterate over is supplied within the script being called.
################################################################################


#### SLURM job control parameters ####

## Options with two comment chars are deactivated

#SBATCH --job-name="parallel-print" #name of the job submitted
#SBATCH --partition=short #name of the queue you are submitting job to
#SBATCH --mail-user=jane.doe@isp.com #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error

module load parallel


#### User-Defined Constants ####

## To specify file to iterate through here, supply its path for iter_file.
## Otherwise set iter_file to a string that is not a valid file name, like "none"
## or "nothing". In this case the file to iterate over should be set in the
## script that is being called.
script="concat_fastqs.sh"
iter_file="nothing"
max_iter=10


#### Executable ####

iter=( $(seq 1 1 ${max_iter}) )

echo
echo "${script}"
echo "Start time:"
date

if [[ -f "$iter_file" ]]; then
    parallel -j $max_iter --delay 1 --joblog parallel_run.log sbatch -t 48:00:00 -N1 -n10 $script $iter_file {} ::: "${iter[@]}"
else
    parallel -j $max_iter --delay 1 --joblog parallel_run.log sbatch -t 48:00:00 -N1 -n10 $script {} ::: "${iter[@]}"
fi

echo
echo "End time:"
date
