#!/bin/bash

## Sequential dispatch script
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
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
## that might list samples, files, genomic regions, etc.
##
## The user can optionally specify the file to iterate over as iter_file. If this
## is set to any string that is not a valid filename, then it assumes that the
## path to the file to iterate over is supplied within the script being called.
################################################################################


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
    for i in "${iter[@]}"; do
        bash "$script" "$iter_file" $i
    done
else
    for i in "${iter[@]}"; do
        bash "$script" $i
    done
fi



echo
echo "End time:"
date
