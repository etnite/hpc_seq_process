#!/bin/bash
set -e

## Clean up Ensembl .gff3 files
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script will get .gff3 files cleaned up and (hopefully) in a valid format.
##
## The input .gff3 file is supplied as a command-line argument. Input can be 
## uncompressed or gzipped. The output file will be uncompressed, with "_cleaned"
## inserted before the extension.
##
## NOTE: The Ensembl .gff3 files do not contain the "chr" prepended to
## chromosome names. They also only contain the IWGSC high-confidence genes,
## in addition to many non-coding RNAs. Here we will just isolate the
## IWGSC high-confidence genes
##
## This script will convert the "Un" unaligned chromosome to "UN", isolate
## IWGSC genes, and sort everything. 
################################################################################


#### SLURM job control #### 

#SBATCH --job-name="call-vars" #name of the job submitted
#SBATCH --partition=short #name of the queue you are submitting job to
#SBATCH --nodes=1 #Number of nodes
  ##SBATCH --ntasks=28  #Number of overall tasks - overrides tasks per node
#SBATCH --ntasks-per-node=22 #number of cores/tasks
#SBATCH --time=36:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=jane.doe@isp.com #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### Executable ####

gff_file=$1

## Remove the ".gff3.gz" for the output name
outname="${gff_file%.*}"
outname=$(echo "$outname" | sed 's/.gff3//')

echo
echo "Writing output to ${outname}_cleaned.gff3"
echo
    
zgrep -v ^"#" "$gff_file" |
    sed -e 's/^Un/UN/' |
    grep "IWGSC" |
    sort -k1,1 -k4,4n > "${outname}_sorted.gff3"

## This code is from the Tabix manual
#(grep ^"#" Triticum_aestivum.IWGSC.41.corrected_UN.gff3; grep -v ^"#" #Triticum_aestivum.IWGSC.41.corrected_UN.gff3 | grep -v "^$" | grep "\t" | #sort -k1,1 -k4,4n) > Triticum_aestivum.IWGSC.41.corrected_UN.sorted.gff3
