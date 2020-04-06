#!/bin/bash

#### SLURM job control #### 

#SBATCH --job-name="untar" #name of the job submitted
#SBATCH --partition=short #name of the queue you are submitting job to
#SBATCH --ntasks=1  #Number of overall tasks - overrides tasks per node
#SBATCH --time=2:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=jane.doe@isp.com #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### User-defined constants ####

tar_archive=""


#### Executable ####

tar -xvf "$tar_archive"
