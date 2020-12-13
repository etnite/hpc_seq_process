#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Split a genome into windows

Brian Ward
brian@brianpward.net
https://github.com/etnite

This script uses pybedtools to split a genome into equally-sized windows. In
addition, it avoids splitting a chromosome inside of genes if the user supplies
a valid .gff gene annotation file. 

The required inputs are:
    1) The samtools-generated fasta index (.fai file) for the reference genome
       fasta file

Optional inputs:
    1) A clean, sorted .gff3 gene annotation file corresponding to the reference
       genome used to create the .fai index file.

The output consists of the following files written to the user-specified out_dir:
    1) A .bed file containing the window intervals
    2) A corresponding regions file as used by samtools/bcftools. This has the format:
         <chromosome>:<1-based start position>-<1-based end position>


NOTES:
    1) pybedtools requires that bedtools is installed in user's PATH
    2) As of this writing, pybedtools does not yet work with python 3.8
"""


#### User-Defined Constants ###################################################


## If there is no .gff3 file, supply any string that doesn't point to a file,
## e.g. "none"
fai_path = '/lustre/project/genolabswheatphg/v2_refseq/iwgsc_refseqv2.0_all_chromosomes.fa.fai'
gff_path = 'none'
out_dir  = '/home/brian.ward/region_files'

## Distance over which to merge together genes that are near each other
## Amount of flanking sequence or "buffer" to add to either end of genic regions
## Desired window size
merge_dist = 5000
flank_dist = 2000
window_size = 100000000


#### Executable ###############################################################


import os
from pathlib import Path
import pybedtools as pbt
import pandas as pd

out_dir = Path(out_dir)
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

## Bedtools makewindows expects a genome file with two tab-delimited columns - 
## The first is chromosome name, and the second is its length. We can get these
## from the first two columns of the .fai index file
## We then create a 1-column dataframe of the chromosome lengths
g_path = out_dir.joinpath('genome_file.txt')
with open(fai_path, 'r') as fai_file, open(g_path, 'w') as g_file:
    for line in fai_file:
        line_split = line.rstrip().split('\t')
        print('\t'.join(line_split[:2]), file = g_file)
genome = pd.read_csv(g_path, sep = '\t', index_col = 0, header = None, )
genome.columns = ['length']

## Set output .bed and regions file paths
bed_path = out_dir.joinpath('regions.bed')
reg_path = out_dir.joinpath('regions.txt')

## If no .gff file is supplied, then we just split up chroms. without regards
## to genic vs. intergenic regions
if not os.path.exists(gff_path):

    ## Split up the genome into windows
    init_wind = pbt.BedTool().makewindows(g = str(g_path), w = window_size).moveto(bed_path)
    
    ## Create corresponding regions file
    with open(reg_path, 'w') as reg_file:
        for line in init_wind:
            start = int(line[1]) + 1
            print(line[0] + ':' + str(start) + "-" + line[2], file = reg_file)
            
        
## If a .gff file is supplied, we tackle the harder problem of only splitting
## chroms between genes
else:
    
    ## Split up the genome into windows
    init_wind = pbt.BedTool().makewindows(g = str(g_path), w = window_size)
        
    ## Now we find all genic regions in the genome
    ## To do this we "flatten" the .gff3 using merge, then add flanking regions
    ## around genic regions
    gff = pbt.BedTool(gff_path)
    genic = gff.merge(d = merge_dist).flank(g = str(g_path), b = flank_dist)
    
    ## Now create a .bed file (and BedTool object) of just the "breaks" between windows
    breaks_path = out_dir.joinpath('breakpoints.bed')
    with open(breaks_path, 'w') as breaks_file:
        for line in init_wind:
                print(line[0] + '\t' + str(int(line[2]) - 1) + '\t' + line[2], file = breaks_file)
    breaks = pbt.BedTool(breaks_path)
            
    ## Now find the closest intergenic region to each breakpoint
    closest = breaks.closest(genic, d = True)
    
    
    '''
    Rather complicated loop to print output .bed and regions .txt files
    We loop through the BedTool object 'closest', and:
       - Start position is set to 0 any time chromosome changes from last iteration
       - If the distance (last col) in 'closest' is 0, we're inside a genic region,
         and we go to closest flank to use as end pos. and next iteration start pos.
       - BUT if we happen to be at end of a chromosome, we always use the position of
         the end of the chromosome, even if it's in a "genic" region (it's probably
         actually inside the user-specified flanking sequence)'''
    chrom = ''
    with open(bed_path, 'w') as bed_file, open(reg_path, 'w') as reg_file:
        for line in closest:
            if line[0] != chrom:
                start_pos = 0
            chrom = line[0]
            if int(line[6]) != 0:
                pos = line[2]
            else:
                borders = [int(line[4]), int(line[5])]
                pos = min(borders, key = lambda x:abs(x - int(line[2])))
            if int(line[2]) == genome.loc[chrom, 'length']:
                pos = line[2]
            start1 = int(start_pos) + 1
            out_list = [chrom, str(start_pos), str(pos)]
            print('\t'.join(out_list), file = bed_file)
            print(chrom + ":" + str(start1) + "-" + str(pos), file = reg_file)
            start_pos = int(pos)
    
    ## Cleanup temporary files
    os.remove(breaks_path)
os.remove(g_path)
    
