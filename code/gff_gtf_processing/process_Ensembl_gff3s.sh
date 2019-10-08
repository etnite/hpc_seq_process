#!/bin/bash
set -e

## Clean up IWGSC .gff3 files
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


gff_file=$1

## Remove the ".gff3.gz" for the output name
outname="${gff_file%.*}"
outname=$(echo "${outname}" | sed 's/.gff3//')

echo
echo "Writing output to ${outname}_cleaned.gff3"
echo
    
zgrep -v ^"#" $gff_file |
    sed -e 's/^Un/UN/' |
    grep "IWGSC" |
    sort -k1,1 -k4,4n > ${outname}_sorted.gff3

## This code is from the Tabix manual
#(grep ^"#" Triticum_aestivum.IWGSC.41.corrected_UN.gff3; grep -v ^"#" #Triticum_aestivum.IWGSC.41.corrected_UN.gff3 | grep -v "^$" | grep "\t" | #sort -k1,1 -k4,4n) > Triticum_aestivum.IWGSC.41.corrected_UN.sorted.gff3
