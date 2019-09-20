# hpc_seq_process

Scripts for processing sequencing data on HPC clusters

## Description

This repository is a collection of scripts for processing sequencing data
(currently just DNA) on high-performance computing (HPC) clusters. It is
currently geared towards clusters using SLURM job schedulers.

## Parallel execution

Many scripts in this repository end with "\_parallel.sh". These scripts are
intended to work on files corresponding to a single sample. They may be used
with a dispatcher script, which coordinates parallel or sequential execution
for multiple samples. These scripts are located in code/parallel_dispatch.
Generally, on clusters arrayer.sh is the preferred script for this purpose.

All of the dispatching scripts simply supply integers within a range starting at
one. This range of integers is then used to subset lines out of a text file,
listing each sample (the path to this file is defined in whatever script is
being controlled).

## Execution on local machines

The code in this repository can be run outside of a HPC setting. In this
use case, parallelizer.sh or looper.sh would be used for parallel (or sequential)
dispatch. Lines within scripts which start with "module" must be commented out.
These lines also show the dependencies for each script. Outside of an HPC
setting, it is usually easiest to use Bioconda environments 
(https://bioconda.github.io/user/install.html) to handle dependencies.
