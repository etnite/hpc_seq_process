#!/bin/bash

## Index a reference genome for alignment with Bowtie2
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script will take a single positional argument - the path to a reference
## genome fasta file, and output a Bowtie2 index with the same base name 
## (i.e. the index file names will have the same format as the input reference 
## file name, without the ".fasta" or ".fa" extension).
##
## NOTE: The Bowtie2 indexer is multi-threaded. Make sure that the number of
## cores specified in the line starting SBATCH -n matches the number of threads
## in the user-defined constants section. Note that on Ceres, the max. number
## of cores per node is 40. Using more cores than this would require MPI, which
## should generally not be necessary in this case.
################################################################################


#### SLURM job control parameters ####

## Options with two comment chars are deactivated

#SBATCH --job-name="bowtie2-index" #name of the job submitted
#SBATCH --partition=mem #name of the queue you are submitting job to
#SBATCH --nodes=1 #Number of nodes
    ##SBATCH --ntasks=10  #Number of overall tasks - overrides tasks per node
#SBATCH --ntasks-per-node=10 #number of cores/tasks
#SBATCH --time=10:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=bpward2@ncsu.edu #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error

module load bowtie2


#### User-defined constants ####

ref_file="/project/genolabswheatphg/v1_refseq/compiled_seqs_noUn/T_aestivum_v1_w_translocs_noUn.fasta"
nthreads=$SLURM_NTASKS


#### Executable ####

echo
echo "Start bowtie2_index_ref.sh"
echo "Start time:"
date

ind_name="${ref_file%.*}"
bowtie2-build --threads $nthreads "${ref_file}" "${ind_name}"

echo
echo "End time:"
date
