#!/bin/bash
set -e

## Calculate depth using samtools depth
##
## This script sums read depths across the genome using samtools depth. It can
## either use a list of bam files supplied in a text file, or else a directory
## containing bam files
##
################################################################################


#### SLURM job control #### 

#SBATCH --job-name="bam-depth" #name of the job submitted
#SBATCH --partition=short #name of the queue you are submitting job to
#SBATCH --nodes=1 #Number of nodes
#SBATCH --ntasks=1  #Number of overall tasks - overrides tasks per node
  ##SBATCH --ntasks-per-node=22 #number of cores/tasks
#SBATCH --time=1:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=jane.doe@isp.com #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### User-Defined Constants ####

## If you want to use all .bam files in a directory, point bam_dir to that
## directory, and set bam_list to a string that doesn't point to any file
bam_list="none"
bam_dir="/project/genolabswheatphg/alignments/ERSGGL_SRW_bw2_bams/SRW_merged_excap_GBS_wholechrom_bw2_bams"

## minimum read quality to be included
## maximum number of reads to count for a position in EACH .bam file
## minimum depth summed ACROSS .bams for a position to be printed
## Maximum gap for bedtools merge
mq_min=20
max_dep=200
min_dep=50
max_gap=100

out_bed="/home/brian.ward/test_out.bed"


#### Executable ####

module load samtools
module load bedtools

## Write the list of .bam files
out_dir=$(dirname "$out_bed")
echo "$out_dir"
if [[ -f "$bam_list" ]]; then
    cp "$bam_list" "$out_dir"/temp_bam_list.txt
else
    ls -1 "$bam_dir"/*.bam > "$out_dir"/temp_bam_list.txt
fi

## Perform depth calculation
samtools depth -Q $mq_min -f "$out_dir"/temp_bam_list.txt |
    awk '{sum=0; for(i=3; i<=NF; i++) { sum+= $i } print($1,$2-1,$2,sum)}' |
    awk -v md=$min_dep 'BEGIN {OFS = "\t"} $4 > md {print($1, $2, $3, $4)}' | head -n 1000 > "$out_dir"/temp.bed
 
bedtools merge -d $max_gap -i "$out_dir"/temp.bed > "$out_bed"

## Run the depth calculation
#samtools depth -Q $mq_min -f "$out_dir"/temp_bam_list.txt | 
#    awk '{for(i=3;i<=NF;i++) sum+=; print($1,$2,sum); sum=0}' | 
#    awk -v md=$min_dep '$3 > md' | awk '$2-lastpos>500 {print lastchr,lastpos,lastdep"\n"$1,$2,$3} {lastchr=$1} {lastpos=$2} {lastdep=$3}' | tail -n +2 | paste - -
