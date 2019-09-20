#!/bin/bash

## Sequential dispatch script
##
## This script is a "universal sequencer" - it's intended function is to
## enable the sequential execution of another script WHEN DIFFERENT PROCESSES
## CAN BE RUN INDEPENDENTLY, such as when running the same process on different
## samples, as is common in bioinformatics.
##
## The script arrayer.sh offers the ability to perform parallel execution in a
## HPC setting, while parallelizer.sh offers the ability to perform parallel
## execution in either HPC or local machines.
##
## This script defines an array of integers, from one to a user-defined maximum.
## Typically, each integer will then be used to subset a single line from a file
## that lists samples (one per line). Note that this file is defined in whatever
## script this one is dispatching (set with the script constant)
################################################################################


#### User-Defined Constants ####

max_iter=10
script="concat_fastqs.sh"


#### Executable ####

iter=( $(seq 1 1 ${max_iter}) )

echo
echo "${script}"
echo "Start time:"
date

for i in "${iter[@]}"; do
    bash "${script}" $i
done

echo
echo "End time:"
date
