#!/usr/bin/env julia

#=
Script for parsing the output of samtools depth. See the argument parser description
below for more details.

Brian Ward
brian@brianpward.net
https://github.com/etnite
=#

using ArgParse
using Statistics

function parse_commandline()

    s = ArgParseSettings(description = "Script for parsing the output of samtools
depth. Reads from stdin and writes to stdout. The user can specify a minimum
mean and/or median read depth to retain each position. The maximum distance
parameter determines the physical distance up to which blocks of contiguous
positions passing the filtering thresholds are merged together, in a manner
similar to bedtools merge. The output is a .bed file of regions passing the
user-specified filters, merged up to the specified distance.\n\n

This script uses the output of samtools depth to 'protect' its input, obviating
the need for some forms of error checking. Samtools depth will also never output a
read depth of more than 8,000 for a position.\n\n

The script may miss one interval, if a window is opened at the end of the last
chromosome, and isn't closed before end-of-file")

    @add_arg_table s begin
        "--min-mean", "-u"
            help = "Float - minimum mean value of reads to include a position in output"
            required = false
            arg_type = Float32
            default = Float32(0)
            range_tester = x -> x >= 0
        "--min-mdn", "-m"
            help = "Float - minimum median value of reads to include a position in output"
            required = false
            arg_type = Float32
            default = Float32(0)
            range_tester = x -> x >= 0
        "--max-dist", "-d"
            help = "Integer - maximum distance over which to merge regions in output .bed file"
            required = false
            arg_type = UInt32
            default = UInt32(0)
            range_tester = x -> x >= 0
    end
    return parse_args(s)
end

function main()

    ## Parse command line arguments
    parsed_args = parse_commandline()
    min_μ = parsed_args["min-mean"]
    min_mdn = parsed_args["min-mdn"]
    dist = parsed_args["max-dist"]

    ## Default mean and median values
    μ = Float32(0)
    mdn = Float32(0)

    ## Initialize position variables
    open_pos = UInt32
    last_pos = Int32(0 - dist)

    ## Iterate through stdin lines
    for line in eachline(stdin)

        ## Get mean and median vals if necessary
        chopped = split(line)
        if min_μ > 0 || min_mdn > 0; depths = parse.(UInt16, chopped[3:end]); end
        if min_μ > 0; μ = mean(depths); end
        if min_mdn > 0; mdn = median(depths); end

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
                    println(stdout, chrom * "\t" * string(open_pos - 1) * "\t" * string(last_pos))
                end ## End check for intial condition
                open_pos = pos
            end ## End check pos - last_pos
            last_pos = pos
        end ## End test mean and median
    end ## End stdin loop
end ## End main() def

main()
