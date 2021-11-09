## Resolve colliding sample names following merging of VCF files
##
###############################################################################


#### User-Defined Constants ###########

in_vcf <- ""
out_vcf <- ""
keep <- "first"
sort_output <- FALSE


#### Executable #######################

line_names <- system(sprintf("bcftools query -l %s", in_vcf), intern = TRUE)

names_df <- data.frame("original" = line_names,
                       "no_prefix" = sub("^[0-9]*:", "", line_names))

dup_names <- unique(no_prefix[duplicated(no_prefix)])
names_df$duplicated <- FALSE
names_df$duplicated[names_df$no_prefix %in% dup_names] <- TRUE

names_df$first <- FALSE
names_df$first[names_df$original == names_df$no_prefix] <- TRUE

original_order <- unique(names_df$no_prefix)

if (keep == "first") {
    names_df <- names_df$first
} else if (keep == "random") {
    singles <- names_df[!names_df$duplicated, ]
    dups <- names_df[names_df$duplicated, ]
    dups <- split(dups, dups$no_prefix)
    for (i in 1:length(dups) {
        shuf <- sample(nrow(dups[[i]]))
        dups[[i]] <- dups[[i]][shuf, ]
        dups[[i]] <- dups[[i]][1, ]
    }
    dups <- do.call("rbind", dups)
    names_df <- rbind(singles, dups)
    names_df <- names_df[match(original_order, names_df$no_prefix), ]
} else {
    stop("Please set keep to either 'first' or 'random')
}


if (sort_output) {
    names_df <- names_df[order(names_df$no_prefix), ]
}
    
write(names_df$no_prefix, file = file.path(dirname(out_vcf, "samp_list.txt")))



system(sprintf("bcftools view %s -S %s -Ob -o %s",
               in_vcf,
               file.path(dirname(out_vcf, "samp_list.txt")),
               out_vcf)
