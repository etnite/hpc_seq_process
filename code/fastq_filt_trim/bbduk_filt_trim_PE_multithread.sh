#!/bin/bash

## Quality filtering and adapter trimming using BBDuk
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This version of the script runs BBDuk using its built-in multi-threading.
## It works well on a single machine/node, but it is difficult to scale beyond
## a single node on an HPC cluster. For that, it is better to use the similar
## parallel version of this script, which runs a single-thread BBDuk process
## on for n different independent sample files. However, BBDuk is quite fast,
## so highly parallel implementation may not always be necessary.
##
## NOTE: bbduk parameters are hard coded in the bbduk.sh call below.
################################################################################


#### SLURM job control ####

## Options with two comment chars are deactivated

#SBATCH --job-name="bbduk-filt" #name of the job submitted
#SBATCH --partition=short #name of the queue you are submitting job to
#SBATCH --nodes=1 #Number of nodes
  ##SBATCH --ntasks=10  #Number of overall tasks - overrides tasks per node
#SBATCH --ntasks-per-node=20 #number of cores/tasks
#SBATCH --time=00:00:30 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=jane.doe@isp.com #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error

module load bbtools


#### User-defined constants ####

fastq_dir="/project/genolabswheatphg/raw_data/wheatCAP_parents"
adapt_fasta="/project/genolabswheatphg/Truseq_paired_adapters.fa"
samples="/project/genolabswheatphg/wheatCAP_samples.tsv"
out_dir="/project/genolabswheatphg/filt_fastqs/wheatCAP_parents"
ncores=$SLURM_NTASKS


#### Executable  ####

echo
echo "Start bbduk_filt_trim_paired_multithread.sh"
echo "Start time:"
date

mkdir "${out_dir}"

## First recursively find all fastq files in fastq_dir
shopt -s globstar nullglob
fastqs=( "$fastq_dir"/**/*.fastq.gz )
#printf '%s\n' "${fastqs[@]}"

## Read the sample names into an array called "samps"
#mapfile -t samps < "${samples}"
samps=( $(cut -f1 "${samples}") )
#printf '%s\n' "${samps[@]}"

## Get length of adapters
ad_len=$(head -n 2 "${adapt_fasta}" | tail -n -1 | wc -c)

for i in "${samps[@]}"; do

    ## Now we have to go through the usual insane shell syntax to get all the unique
    ## lanes for sample i
    lanes=($(printf '%s\n' "${fastqs[@]}" | grep "$i" | sed 's/_R[12].*$//' | sort -u))
    #echo "***"    
    #printf '%s\n' "${lanes[@]}"

    ## Loop through each lane
    for j in "${lanes[@]}"; do

	base=$(basename "${j}")
	#echo "samp: ${i}"
        echo "lane: ${j}"

        ## Assign yet another new array, holding each of the two mated fastqs for lane i    
        mates=($(printf '%s\n' "${fastqs[@]}" | grep $j))    
        fq1="${mates[0]}"
        fq2="${mates[1]}"
        #echo "***"
        #echo $fq1
        #echo $fq2

        ## Run BBDuk
        ## Setting maq will remove entire reads based on their average qual
        ## maq=10 average error is 10%
        ## maq=13 ~5%
        ## maq=15 ~3%
        ## maq=20 1%
        bbduk.sh -Xmx10g in1="${fq1}" in2="${fq2}" \
            out1="${out_dir}"/"${base}"_R1.fastq.gz out2="${out_dir}"/"${base}"_R2.fastq.gz \
            ref="${adapt_fasta}" ktrim=r k=$ad_len mink=10 hdist=3 hdist2=1 ftm=5 maq=13 minlen=75 tpe tbo
    done
done

echo
echo "End time:"
date
