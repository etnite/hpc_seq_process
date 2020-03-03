#!/usr/bin/env julia

#=
Script for parsing the output of samtools depth. See the argument parser description
below for more details.
=#

using ArgParse
using Statistics

function parse_commandline()

    s = ArgParseSettings(description = "Script for parsing the output of samtools depth. The user can specify a minimum
mean and/or median read depth to retain each position. The maximum distance
parameter determines the physical distance up to which blocks of contiguous
positions passing the filtering thresholds are merged together, in a manner
similar to bedtools merge. The output is a .bed file of regions passing the
user-specified filters, merged up to the specified distance.\n\n

This script uses the output of samtools depth to 'protect' its input, obviating
the need for some forms of error checking\n\n
    
The script may miss one interval, if a window is opened at the end of the last
chromosome, and isn't closed before end-of-file")

    @add_arg_table s begin
        "--min-mean", "-u"
            help = "Float - minimum mean value of reads to include a position in output"
            required = true
            arg_type = Float32
        "--min-mdn", "-m"
            help = "Float - minimum median value of reads to include a position in output"
            required = true
            arg_type = Float32
        "--max-dist", "-d"
            help = "Integer - maximum distance over which to merge regions in output .bed file"
            required = true
            arg_type = UInt32
    end
    return parse_args(s)
end

function main()

    ## Parse command line arguments
    parsed_args = parse_commandline()
    min_μ = parsed_args["min-mean"]
    min_mdn = parsed_args["min-mdn"]
    dist = parsed_args["max-dist"]

    ## Initialize position variables
    open_pos = UInt32
    last_pos = Int32(0 - dist)

    ## Open std input stream; iterate through lines
    for line in eachline(stdin)

        ## Get mean and median vals
        chopped = split(line)
        depths = parse.(UInt16, chopped[3:end])
        μ = mean(depths)
        mdn = median(depths)

        ## Test if the position's mean and median depths are above thresholds
        if μ >= min_μ && mdn >= min_mdn
            chrom = chopped[1]
            pos = parse(UInt32, chopped[2])

            ## Iteratively update the position opening a window (open_pos)
            ## If starting, open_pos is set to first position satisfying mean/median filters
            ## Otherwise, if the current position is > dist away from last_pos
            ## we write out the previous open_pos and the last_pos to the .bed
            ## file, and then set open_pos to current position
            if abs(pos - last_pos) > dist
                if last_pos >= 0
                    println(chrom * "\t" * string(open_pos - 1) * "\t" * string(last_pos))
                end ## End check for intial condition
                open_pos = pos
            end ## End check pos - last_pos
            last_pos = pos
        end ## End test mean and median
    end ## End stdin loop
end ## End main() def

main()
