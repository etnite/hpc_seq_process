#!/bin/bash

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

bam_list_file=""
max_dep=100
out_file=""


#### Executable ####

module load samtools

samtools depth -m $max_dep -f "$bam_list_file" |
	head -n 100 |
	./samtools_depth_summstats.py |
	gzip -c > "$out_file"
