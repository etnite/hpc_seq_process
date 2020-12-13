#!/usr/bin/env julia

#=

=#

using CSV
using DataFrames
using BGZFStreams
using Statistics

const bed_file = "/home/brian/Downloads/GDrive_downloads/chr_translate.bed"
const in_vcf_file = "/home/brian/Downloads/GDrive_downloads/test.vcf.gz"
const out_vcf_file = "/home/brian/Downloads/GDrive_downloads/test_out.vcf.gz"

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


for line in eachline(istream)
    if startswith(line, "##contig")
        continue
    elseif startswith(line, "#CHROM")


end
##contig=<ID=20,length=62435964
