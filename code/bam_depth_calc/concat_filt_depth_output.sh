#!/bin/bash

## Merge samtools_depth_parallel.sh output
##
## This is just a simple script to turn the output of the associated script
## samtools_depth_parallel.sh into a single, more compact .bed file. The bedtools
## merge tools supports adding in various notations in an additional fourth
## column. This gets a little tricky to interpret though. In the example below,
## we filter our samtools depth output by the median read depth at each position.
## Then we merge together everything separated by less than the length of a read
## (in this case 150bp). In the output .bed file we also calculate the mean of
## column 4, which is the number of reads at that position. Note that this IS NOT
## the mean read depth at a particular position - it is the mean of the read
## depths across the entire interval listed in the output .bed file.
################################################################################


#### User-Defined Constants ####

in_dir=""
out_bed=""


#### Executable ####

cat "${out_dir}/!(*header).txt" |
    awk 'BEGIN{OFS="\t"} $7 >= 4 {print $1, $2, $3}' |
    bedtools merge -i stdin -d 150 -c 4 -o mean > "$out_bed"
