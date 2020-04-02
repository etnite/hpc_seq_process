#!/bin/bash

## Test script to print out contents of a file using parallelizer.sh

date
echo "input: $1" 
echo "number of nodes: ${SLURM_NNODES}"
echo "nodes: ${SLURM_JOB_NODELIST}"
echo "number of cores: ${SLURM_NPROCS}"
echo "number of tasks: ${SLURM_NTASKS}"
echo "job ID: ${SLURM_JOBID}"
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
