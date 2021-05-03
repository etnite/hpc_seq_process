#!/usr/bin/env julia

#=
Convert VCF file from chromosome parts or sub-feature coordinates to full
chromosome coordinates

This script takes as input a VCF file (compressed with gzip/bgzip or uncompressed)
and updates the coordinates based on a .bed file. It is intended primarily to
convert a VCF file with coordinates in terms of chromosome pieces, scaffolds,
contigs, etc. into full chromosome pseudomolecule coordinates.

This is useful because some programs, notably GATK, cannot handle positions
higher than approximately 530Mb.

The .bed file directs the coordinate conversion. This is a tab-delimited file
(no header) with the following columns:
  1) Chromosome name (feature will be renamed to this in output)
  2) Feature start position (0-based)
  3) Feature end position (1-based)
  4) Feature name (current name of feature in VCF)

Example of first few rows of a .bed file (using "|" to denote tabs):

1A |         0 | 471304005 | chr1A_part1
1A | 471304005 | 594102056 | chr1A_part2
1B |         0 | 438720154 | chr1B_part1

=#

using CSV
using DataFrames
using BGZFStreams
using Statistics

bed_file = "/home/brian/Downloads/GDrive_downloads/chr_translate.bed"
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
maxes = groupby(bed, :Column1)
maxes = combine(maxes, :Column3 => maximum => :max)

## Create vector of new contig IDs
contigs = []
for i in 1:nrow(maxes)
    push!(contigs, "##contig=<ID=" * string(maxes.Column1[i]) * ",length=" * string(maxes.max[i]) * ">")
end

## Create dictionary to translate between sequences
seq_dict = Dict{AbstractString, Tuple{AbstractString, UInt64}}()
for i in 1:nrow(bed)
    seq_dict[bed.Column4[i]] = (bed.Column1[i], bed.Column2[i])
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
        orig_chrom = line[1]
        line[1] = seq_dict[orig_chrom][1]
        line[2] = string(seq_dict[orig_chrom][2] + parse(UInt64, line[2]))
        #new_chr = seq_dict[line[1]][1]
        #new_pos = seq_dict[line[1]][2] + parse(UInt64, line[2])
        #line = [new_chr; string(new_pos); line[3:length(line)]]
        println(ostream, join(line, "\t"))
    end
end

println("Finished!")
close(istream)
close(ostream)
