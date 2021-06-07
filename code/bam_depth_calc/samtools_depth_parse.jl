#!/usr/bin/env julia

#=
Script for parsing the output of samtools depth.

This script reads from standard in and writes to standard out.

The format of samtools depth is a tab-delimited file with the following columns:
   1) Chromosome
   2) Position
   3-end) All remaining columns as read depths of each genotype at the specified position

The output of this script is a tab-delimited file with the following columns:
   1) Chromosome
   2) Position - 1 (making this sort of a quasi-bed file format)
   3) Position
   4) Sum of read depths
   5) Minimum of read depths
   6) Maximum of read depths
   7) Median of read depths
   8) Mean of read depths
   9) Standard deviation of read depths

Brian Ward
brian@brianpward.net
https://github.com/etnite
=#

using Statistics


## Iterate through stdin lines
for line in eachline(stdin)
  
    chopped = split(line)
    chrom = chopped[1]
    pos = parse(UInt32, chopped[2])
    depths = parse.(UInt16, chopped[3:end])

    ## Print to stdout
    println(stdout, join(
        [chrom,
        string(pos - 1),
        string(pos),
        string(sum(depths)),
        string(minimum(depths)),
        string(maximum(depths)),
        string(median(depths)),
        string(mean(depths)),
        string(std(depths))],
        "\t"
    ))
    flush(stdout)

end
