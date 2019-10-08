#!/bin/bash
set -e

## Convert a .gff3 file to a flattened .bed file
##
## This script will take a .gff3 file as input, sort it, isolate only lines that
## represent genes (i.e. get rid of subfeatures), and then will "flatten" it.
## The input .gff3 file can be gzipped or uncompressed.
##
## In this context, flattening means that overlapping genes will be merged into
## single intervals. The output is a .bed file representing these flattened
## intervals. Bedtools is required to perform the flattening ("merging" in
## bedtools documentation)
##
## NOTES: I don't believe that bedtools merge can read from standard input, so
## a temporary file is created and deleted. Also, note that .gff3 files
## have feature start positions in 1-based coordinates, while these are converted
## to 0-based coordinates in the output .bed file. Bedtools seems to recognize
## gff3 files even if they don't have the .gff3 extension
################################################################################

gff_file=$1

## Remove the ".gff3.gz" for the output name
outname="${gff_file%.*}"
outname=$(echo "${outname}" | sed 's/.gff3//')
outdir=$(dirname "${outname}")

tmpfile=$(mktemp -p "${outdir}")

## Get only genes from .gff3 and sort 
zgrep -P "\tgene\t" "${gff_file}" |
    sort -k1,1 -k4,4n -k5,5n > "${tmpfile}"

bedtools merge -i "${tmpfile}" > "${outname}"_flat.bed

rm "${tmpfile}"
