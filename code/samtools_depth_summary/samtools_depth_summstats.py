#!/usr/bin/env python3

'''
Calculate summary stats for samtools depth output

This script parses the output of the command:

    samtools depth [input_bam|-f bam_filelist]

The output of this command is a tab-delimited file with first column
chromosome, second column position, and all subsequent columns consisting of
read depths for individual SAM/BAM files at that position.

The script takes one positional argument - the minimum mean read depth. If
mean read depth is below this number, the position is not output

The script writes to stdout, so the output can then be piped to a file, or
gzipped if desired.
'''

import fileinput
import sys
import numpy as np

## Get min_mean val
#if len(sys.argv) > 1:
#    min_mean = int(sys.argv[1])
#else:
#    min_mean = 0
min_mean = 5

## Read input line from stdin, split SNP info from read depths
line_count = 0
for line in fileinput.input():
    
    ## Split up the line
    line_split = line.rstrip().split('\t')
    pos_info = line_split[0:2]
    depths = np.array(line_split[2:]).astype(np.int)

    ## Calculate summary stats
    mean_dep = np.mean(depths)
    if mean_dep >= min_mean:
        pos_summ = [np.sum(depths), mean_dep, np.median(depths), np.max(depths)]

        ## Write output
        out_list = pos_info + pos_summ
        print(*out_list, sep = '\t', end = '\n')

    ## Print progress
    line_count += 1
    if line_count % 10000 == 0:
        print(str(line_count/1000) + 'K positions processed', file = sys.stderr)
