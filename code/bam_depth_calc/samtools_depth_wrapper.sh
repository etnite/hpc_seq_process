#!/bin/bash

## Wrapper script for samtools_depth_parse.jl
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script handles input and output streaming for samtools_depth_parse.jl
## and allows for running the parser in parallel across a list of genome regions
##
## The samtools depth command outputs a large matrix, where: 
##   * each cell is a tally of read depths
##   * each row is a chromosome/position combination
##   * each column is a sample
## samtools_depth_parse.pl parses this output, to discard positions where 
## mean and/or median read depth falls below user-specified thresholds. It then
## merges contiguous positions passing the thresholds into regions. A maximum
## distance over which to merge positions into regions can be specified. Finally,
## following merging any regions that are smaller than a user-defined value are
## discarded.
##
## Regions of any size can be used (using for instance bedtools makewindows).
## But keep this in mind:
##   1) regions that are smaller than whole chroms can cause artifacts at
##      boundaries (e.g. intervals being deleted for being too small)
##   2) A file defining regions must be in samtools format chr:from_bp-to_bp, or
##      optionally just chromosome, one region per line. A .bed file defining 
##      regions won't work here, as GNU parallel must take its input source 
##      line-by-line from a file
################################################################################


#### SLURM job control #### 

#SBATCH --job-name="par-dep" #name of the job submitted
#SBATCH --partition=short #name of the queue you are submitting job to
#SBATCH --nodes=1 #Number of nodes
    ##SBATCH --ntasks=1  #Number of overall tasks - overrides tasks per node
#SBATCH --ntasks-per-node=22 #number of cores/tasks
#SBATCH --time=08:00:00 #time allocated for this job hours:mins:seconds
#SBATCH --mail-user=jane.doe@isp.com #enter your email address to receive emails
#SBATCH --mail-type=BEGIN,END,FAIL #will receive an email when job starts, ends or fails
#SBATCH --output="stdout.%j.%N" # standard out %j adds job number to outputfile name and %N adds the node name
#SBATCH --error="stderr.%j.%N" #optional but it prints our standard error


#### User-Defined Constants ####

## Can provide a text file listing bam files, one per line
## Or supply a string that isn't a filename for bam_list_file
## And supply a directory to look for .bam files
bam_list_file="none"
bam_dir="/project/genolabswheatphg/alignments/ERSGGL_SRW_bw2_bams/SRW_merged_excap_GBS_wholechrom_bw2_bams_mq20_filt"

## Path to reference genome .fasta file (for finding chromosome names)
ref_gen="/project/genolabswheatphg/v1_refseq/whole_chroms/Triticum_aestivum.IWGSC.dna.toplevel.fa"

## The following parameters are for:
## 1) Maximum read depth PER SAMPLE, after which samtools depth stops counting
## 2) Minimum mapping quality to consider a read
## 3) Min mean depth to include a position, calculated ACROSS SAMPLES
## 4) Min median depth to include a position, calculated ACROSS SAMPLES 
## 5) Maximum distance between features in output .bed file. If two features are closer together
##    than this number, they are merged into a single feature (similar to -d parameter of
##    bedtools merge)
## 6) Minimum feature size in output .bed file - features smaller than this are discarded
max_dep=200
min_qual=20
min_mean=0
min_median=5
max_dist=500
min_size=50

## Path to output .bed file
out_bed="/project/genolabswheatphg/SRW_test_phg/find_ref_ranges/SRW_mdn5_dist500_size50_dp.bed"


#### Executable ####

module load samtools
module load parallel

echo
echo "Start time:"
date

## Create the list of .bam files
out_dir=$(dirname "$out_bed")
mkdir -p "$out_dir"
mkdir "$out_dir"/region_depths
if [[ -f "$bam_list_file" ]]; then
    cp "$bam_list_file" "$out_dir"/temp_bam_list.txt
else
    realpath "$bam_dir"/*.bam > "$out_dir"/temp_bam_list.txt
fi

## One method of making a regions file (in this case of whole chroms)
if [[ ! -f "${ref_gen}".fai ]]; then samtools faidx "${ref_gen}"; fi
cut -f 1 "${ref_gen}".fai > "$out_dir"/regions.txt

## Run depth calculation; Lots going on here
## Sometimes easier to give parallel a function, but in this case there would
## be many positional args to supply
parallel -j $SLURM_NTASKS -a "$out_dir"/regions.txt \
    "samtools depth -r {} -m $max_dep -Q $min_qual -f ${out_dir}/temp_bam_list.txt |
    ./samtools_depth_parse.jl -u $min_mean -m $min_median -d $max_dist |
    awk '(\$3 - \$2) >= $min_size' > ${out_dir}/region_depths/{}.bed" #::: "${chroms[@]}"

#    ./samtools_depth_summstats.py |
#    gzip -c > "$out_bed"

## Combine region .bed files and delete temporary BAM list file
cat "$out_dir"/region_depths/*.bed | sort -k 1,1 -k 2,2n > "$out_bed"
#rm "$out_dir"/temp_bam_list.txt

echo
echo "End time:"
date
