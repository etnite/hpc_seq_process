#!/bin/bash
set -e
#source /home/gbg_lab_admin/miniconda3/bin/activate bwa_align_call


## Call variants in parallel
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script performs variant calling using BCFTools, running the bcftools
## mpileup and call commands in parallel on a list of regions. The user can supply
## a list, which should NOT be .bed file. Instead, it should be in samtool's
## usual regions format:
##
##    chromosome:start_pos-end_pos
##
## One region per line. If no region file is supplied, then the script creates
## one from the reference genome, running each chromosome separately.
##
## A minimum mapping-quality threshold may be set. For bowtie2, probability of
## an incorrect alignment (p) is related to the mapping quality value (Q) by:
##
##   p = 10 ^ -(Q/10)
##
## A mapq score of 5 translates to a ~32% chance of misalignment
## A mapq score of 10 translates to a 10% chance
## A mapq score of 20 translates to a 1% chance, etc.
##
## Currently the script only uses reads mapped in "proper pairs". This can be disabled
## by deleting the "--rf 2" in the bcftools mpileup call.
##
## NOTE: Can save the intermediate, mpileup-generated BCFs for each chrom by setting
## save_pile_out to "true". This is slower, as we then must compress each
## of the mpileup BCFs, and then uncompress to feed into bcftools call. Saving uncompressed
## mpileup BCFs is probably not a great idea, because they will be very large.
##
## The script outputs a BCF file rather than a VCF, as BCF files are faster to
## process downstream using bcftools.
##
## The script creates a directory for indivual chromosome variant-only BCFs, but this is removed
## at the end, as the final BCF should contain all the same information.
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


#### User-defined constants ####

bams_list="/project/genolabswheatphg/alignments/ERSGGL_SRW_bw2_bams/SRW_merged_excap_GBS_wholechrom_bw2_bams"
ref_gen="/project/genolabswheatphg/v1_refseq/whole_chroms/Triticum_aestivum.IWGSC.dna.toplevel.fa"
regions_file=""
out_bcf="/project/genolabswheatphg/variants/SRW/ERSGGL_SRW_merged_excap_GBS_wholechr_bw2.vcf.gz"
ncores=22
mq_val=20
save_pile_out="false"


#### Executable ####

module load parallel
module load samtools
module load bcftools

echo
echo "Start time:"
date

## Grab first letter of save_pile_out
save_pile_out="${save_pile_out:0:1}"

## Create output directory, cd to it, create subdirectories for single-chrom BCFs
out_dir=$(dirname "$out_bcf")
mkdir -p "$out_dir"
cd "$out_dir"
mkdir chrom_var_bcfs

## Sanity check on save_pile_out
if [[ "$save_pile_out" == [Tt] ]]; then
    mkdir mpi_bcfs
elif [[ "$save_pile_out" == [Ff] ]]; then
    echo "User elected not to save mpileup BCF files"
else
    echo "Please supply 'true' or 'false' for save_pile_out"
    exit 1;
fi

## Store list of regions in array
## If a file of regions is supplied we use it; otherwise we use a list of the
## different chromosomes generated from the reference genome
if [[ ! -f "$regions_file" ]]; then
    regions=( $(cat "$regions_file") )
else
    if [[ ! -f "${ref_gen}.fai" ]]; then 
        samtools faidx "$ref_gen"
        regions=( $(cut -f 1 "${ref_gen}.fai") )
    fi
fi

## Then run variant calling pipeline on chromosomes in parallel
## There are a lot of options for mpileup and call - using some defaults for now
## By default, mpileup will skip reads that are unmapped, secondary,
##   PCR duplicates, or that failed QC.
if [[ "$save_pile_out" == [Tt] ]]; then
    time parallel -j $ncores bcftools mpileup -Ob -q $mq_val --rf 2 -a FORMAT/DP,FORMAT/AD -f "$ref_gen" -b "$bams_list" -r {} -o mpi_bcfs/chrom_{}.bcf ::: "${regions[@]}"
    time parallel -j $ncores bcftools call -mv -Ou mpi_bcfs/chrom_{}.bcf -o chrom_var_bcfs/chrom_{}.bcf ::: "${regions[@]}"
else
    time parallel -j $ncores "bcftools mpileup -Ou -S $samples -q $mq_val --rf 2 -a FORMAT/DP,FORMAT/AD -f $ref_gen -b $bams_list -r {} | bcftools call -mv -Ou -o chrom_var_bcfs/chrom_{}.bcf" ::: "${regions[@]}"
fi

## Create list of single-chromosome variant BCFs
printf '%s\n' chrom_var_bcfs/*.bcf > chrom_var_bcfs/bcf_list.txt

## Concatenate all chromosome BCFs into single gzipped VCF
## Construct SNP IDs
## Normalize indels; discard all but one overlapping SNP, all but one overlapping indel
bcftools concat --no-version --threads $ncores --file-list chrom_var_bcfs/bcf_list.txt -Ou |
    bcftools annotate --set-id +'S%CHROM\_%POS' -Ou |
    bcftools norm -f ${ref_gen} -Ob > "$out_bcf"
 
bcftools index -c "$out_bcf"

rm -rf "$chrom_var_bcfs"

echo
echo "End time:"
date
