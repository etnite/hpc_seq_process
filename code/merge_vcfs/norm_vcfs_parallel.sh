#!/bin/bash

## Normalize multiple VCF/BCF files
##
## Brian Ward
## brian@brianpward.net
## https://github.com/etnite
##
## This script will normalize a list of VCF or BCF files. It requires that bcftools
## be installed in the user's PATH.
##
## Normalization may be required prior to merging together VCFs, depending upon
## the SNP calling method that was used. For instance, Tassel-GBS sets the major
## allele as reference, and therefore reference alleles won't match between
## different VCF files created using Tassel.
##
## The script takes two positional inputs - the first is a list of VCF files to
## merge. The second is an integer specifying which line of the file to process.
## If 0 is supplied, then all listed VCFs are normalized sequentially.
##
## In addition to the normalized BCF files, the script will also output a list
## containing the names of the normalized files in the output directory.
##
## NOTES ON NORMALIZATION
##   The form of normalization used here will switch alleles so that the reference
##   allele matches the actual allele at that position in the reference genome.
##   This must be exactly the same reference genome used for the SNP calling in
##   the first place. If neither allele at a position matches the actual REF allele,
##   then the REF allele is created as an extra "phantom" allele. Only the GT and
##   AC allele count fields in the VCF are updated - other fields remain unchanged.
################################################################################


##### User-Defined Constants ##########

## Path to place normalized BCF files
out_dir="/home/gbg_lab_admin/Array_60TB/Wheat_GBS/all_nurseries_Nov_2021/raw_VCF/nurseries_norgrains_merge"

## Path to reference genome
ref_genome="/home/gbg_lab_admin/Array_60TB/GBS_Reference_Genomes/Ensembl_v41_IWGSC_v1.0/Triticum_aestivum.IWGSC.dna.toplevel.fa"

## true/false specifying whether to only retain biallelic SNPs following normalization
biallelic_snps_only="true"


##### Executable #######################

module load singularity/3.7.1
module load bcftools

echo
echo "Start time:"
date

## Get command line args
in_list="$1"
array_ind="$2"

## Get first letter of true/false options
biallelic_snps_only="${biallelic_snps_only:0:1}"

## Set up the output directories
mkdir -p "$out_dir"
if [[ $array_ind -eq 0 || $array_ind -eq 1 ]]; then
    rm -f "${out_dir}/normed_list.txt"
else
    sleep 20s
fi

## Construct an array - either the entire input list, or else one line of it
if [[ $array_ind -eq 0 ]]; then
    mapfile -t vcf_array < "$in_list"
else
    n_toss=$(( $array_ind - 1 ))
    mapfile -t -n $array_ind -s $n_toss vcf_array < "$in_list"   ## Should work for -s 0
fi

## Now loop through our array of input file(s)
for i in "${vcf_array[@]}"; do
    if [[ ! -f "${i}.csi" ]]; then bcftools index -c "$i"; fi

    ## Get the basename of the input VCF/BCF
    ## Then we create a corresponding filename for after normalization
    ## This will always be a .bcf file and may (possibly) need its extension changed
    base=$(basename "$i")
    normed_file="${out_dir}/${base}"
    normed_file="${normed_file%.*}"
    normed_file="${normed_file%.*}"
    if [[ "$biallelic_snps_only" == [Tt] ]]; then
        normed_file="${normed_file}_normed_biallelic_snps.bcf"
    else
        normed_file="${normed_file}_normed.bcf"
    fi
    echo "$normed_file" >> "${out_dir}/normed_list.txt"

    echo
    echo "Normalizing file $i"

    if [[ "$biallelic_snps_only" == [Tt] ]]; then
        bcftools norm "$i" \
            --fasta-ref "$ref_genome" \
            --check-ref s \
            --multiallelics +both \
            --output-type u |
            bcftools view - -m2 -M2 -v snps \
                --output-type b \
                --output "$normed_file"
    else
        bcftools norm "$i" \
            --fasta-ref "$ref_genome" \
            --check-ref s \
            --multiallelics +both \
            --output-type b \
            --output "$normed_file"
    fi

    echo
    echo "Indexing file $normed_file"
    bcftools index -c "$normed_file"
done

sort "${out_dir}/normed_list.txt" -o "${out_dir}/normed_list.txt"


echo
echo "End time:"
date
