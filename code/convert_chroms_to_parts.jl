#!/usr/bin/env julia

#=
Convert VCF file from chromosome pseudomolecule coordinates to feature coordinates

This script takes as input a VCF file (compressed with gzip/bgzip or uncompressed)
and updates the coordinates based on a .bed file. It is intended primarily to
convert a VCF file with coordinates in terms of full chromosome pseudomolecules
into smaller scaffolds, contigs, parts, etc.

The .bed file directs the coordinate conversion. This is a tab-delimited file
(no header) with the following columns:
  1) Chromosome name (current name of chromosome in VCF)
  2) Feature start position (0-based)
  3) Feature end position (1-based)
  4) Feature name (name will be updated to this in output VCF)

Example of first few rows of a .bed file (using "|" to denote tabs):

1A |         0 | 471304005 | chr1A_part1
1A | 471304005 | 594102056 | chr1A_part2
1B |         0 | 438720154 | chr1B_part1

NOTE: features cannot be overlapping. If this occurs, any variants occuring within
multiple features will be discarded with a warning. Similarly, variants that do
not occur in any feature are also discarded with a warning.
=#

using CSV
using DataFrames
using BGZFStreams
using Statistics

bed_file = "/home/brian/Downloads/GDrive_downloads/test_chr_convert.bed"
in_vcf_file = "/home/brian/Downloads/GDrive_downloads/test.vcf.gz"
out_vcf_file = "/home/brian/Downloads/GDrive_downloads/test_out.vcf"

## Here we open input and output streams
## If filename ends with "gz" we read and/or write using bgzip
if endswith(in_vcf_file, "gz")
    istream = BGZFStream(in_vcf_file, "r")
else
    istream = open(in_vcf_file, "r")
end

if endswith(out_vcf_file, "gz")
    ostream = BGZFStream(out_vcf_file, "w")
else
    ostream = open(out_vcf_file, "w")
end

## Read in the .bed file and find max. position of each chromosome
bed = CSV.read(bed_file, DataFrame, header = false, delim = "\t")
rename!(bed, ["chrom", "start", "end", "feature"])

## Create vector of new contig IDs
contigs = []
for i in 1:nrow(bed)
    push!(contigs, "##contig=<ID=" * string(bed.feature[i]) * ",length=" * string(bed.end[i] - bed.start[i]) * ">")
end

#### Main Loop ####
## Remove existing contig lines; dump new contig lines before #CHROM line
## Update chromosome names and positions in data lines
println("Converting VCF position coordinates...")
for line in eachline(istream)
    if startswith(line, "#")
        if startswith(line, "##contig")
            continue
        elseif startswith(line, "#CHROM")
            for contig in contigs; println(ostream, contig); end
            println(ostream, line)
        else
            println(ostream, line)
        end
    else
        line = split(line, "\t")
        pos = parse(UInt64, line[2])

        bed_row = bed[(bed.chrom .== line[1]) .& (bed.start .< pos) .& (bed.end .>= pos), :]
        if nrow(bed_row) == 0
            println("WARNING: Variant at chrom: " * line[1] * " pos: " * line[2] * " does not occur in any listed feature - skipping")
        elseif nrow(bed_row) > 1
            println("WARNING: Variant at chrom: " * line[1] * " pos: " * line[2] * " occurs in multiple features - skipping")
        else
            line[1] = bed_row[1, :feature]
            line[2] = string(pos - bed_row[1, :start])
            println(ostream, join(line, "\t"))
        end
    end
end

println("Finished!")
close(istream)
close(ostream)
