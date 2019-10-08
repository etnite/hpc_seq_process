#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Joseph Schroedl
#jeschroe@ncsu.edu
#https://github.com/corndog2000

## See argument parser epilog for program description

import os
import argparse
from argparse import RawTextHelpFormatter
import csv

def parse_arguments():

    ## Parse command line arguments
    parser = argparse.ArgumentParser(prog = "MakeRegionsFile.py",
    description="Create a file with number ranges based on input.", 
    formatter_class=RawTextHelpFormatter, 
    epilog='''
    This program takes as input an index of a fasta file created with
    samtools faidx (typically a pseudomolecule-level reference genome). It will 
    then "fragment" this fasta file, producing an output BED file with each
    fasta file entry divided into intervals of a size of the user's choosing. 
    The output file is in BED format so the interval starting number is 
    zero-based and the ending number is one-based.

    The output .bed file will be located in "output.bed", in the same directory
    as the input .fai file.
       
    Example output:
        1A	0	10
        1A	10	20
        1A	20	30
    ''')

    parser.add_argument("input_file",
        help = "Path to the input .fai file", type = str)
    parser.add_argument("interval", 
        help = "Number that specifies the size of each range of numbers")
    return parser.parse_args()

def main():

    ## Assign command line args to global variables
    args = parse_arguments()
    input_file = args.input_file
    interval = int(float(args.interval))
    
    ## Create the output file in the same directory as the input file
    output_file = os.path.dirname(input_file)
    
    if input_file[len(input_file) - 1] != "/" and len(os.path.dirname(input_file)) > 0:
        output_file = output_file + "/"
    
    output_file = output_file + "output.bed"

    #print(output_file)


    ## Check for all required arguments
    kill_prog = False

    if input_file is None:
        print("\033[91m" + "Error: Missing input_file argument. Please follow usage instructions." + "\x1b[0m")
        kill_prog = True
    elif os.path.exists(input_file) is False:
        print("\033[91m" + "Error: input_file does not exist" + "\x1b[0m")
        kill_prog = True
    if interval is None:
        print("\033[91m" + "Error: Missing interval argument. Please follow usage instructions." + "\x1b[0m")
        kill_prog = True

    ## If any of the requirements are not met then end the program
    if kill_prog is True:
        exit()

    with open(input_file, "r", newline = "") as csvfile:
        reader = csv.reader(csvfile, delimiter = "\t")

        with open(output_file, "w+") as output:

            ## Get the chromosome and number of base pairs from the input file
            for row in reader:
                chrom = row[0]
                total_pairs = int(row[1])
                written_pairs = 0

                ## Write all the rows based on the interval
                while written_pairs < total_pairs:
                    ## Add the chromosome to the eventual output
                    to_write = f"{chrom}\t"
                    ## Add the current starting value of the interval to the output
                    
                    ## Uncomment the line below to make the interval start 1-based
                    #to_write = to_write + f"{written_pairs + 1}\t"
                    
                    ## Comment out the line below to make the interval start 1-based
                    to_write = to_write + f"{written_pairs}\t"
                    
                    ## If the remaining base pairs don't fit into a whole new row then add them onto the current last one
                    if ((total_pairs - written_pairs) / interval) < 1:
                        written_pairs = written_pairs + (total_pairs - written_pairs)
                    else:
                        written_pairs = written_pairs + interval

                    ## Write the interval end this is 0-based    
                    to_write = to_write + f"{written_pairs}\n"

                    ## Write the new line to output file
                    output.writelines(to_write)

    print("Done")


if __name__ == "__main__":
    main()
