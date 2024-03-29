## Calculate and plot stats from a BCF/VCF file
##
## This script takes as input one BCF or VCF file, and will then output both
## variant- and sample-wise summary statistics, as well some histograms.
##
## BCFTools is required to extract information from the input file. Therefore the
## program must either be installed system-wide or else R/Rstudio must be started
## from within a Conda environment containing BCFTools.
##
## NOTE that this script only uses biallelic variants - the summation of heterozygous
## calls isn't straightforward for variants with > 2 alleles.
##
## I have tried this script out on Linux and MacOS - no idea if it would work on
## Windows.
##
## Finally, the steps in which data is pulled out of the input file can be slow
## for large files. The process will be significantly faster for BCF files compared
## to VCF files.
################################################################################
rm(list = ls())
library(R.utils)
library(data.table)
library(psych)
library(ggplot2)


#### User-Defined Constants ####################################################

vcf_file <- "/autofs/bioinformatics-ward/VCF_generation/OSU_samps_for_Nelly_Apr2022/OSU_samples_for_Nelly_April2022_raw.bcf"
wkdir <- "/autofs/bioinformatics-ward/VCF_generation/OSU_samps_for_Nelly_Apr2022/Nelly_raw_BCF_stats+plots"

## Remove BCFTools-generated files at the end (T/F)?
delete_interfiles <- TRUE


#### Executable ################################################################

dir.create(wkdir, recursive = TRUE)
setwd(wkdir)
dir.create("plots")


#### Generate intermediate files ####

## Generate the relevant variants info .txt file
system(sprintf("bcftools view %s --min-alleles 2 --max-alleles 2 -Ou | 
                bcftools +fill-tags -Ou - -- -t MAF,F_MISSING,AC_Het,NS |
                bcftools +fill-tags -Ou - -- -t 'DPsum:1=int(sum(FORMAT/DP))' |
                bcftools query -f '%%CHROM\t%%POS\t%%INFO/DPsum\t%%INFO/MAF\t%%INFO/F_MISSING\t%%INFO/AC_Het\t%%INFO/NS\n' |
                gzip -c > %s", vcf_file, "var_stats.txt.gz"))

## And the samples info .txt file
system(sprintf("bcftools stats -s - %s | 
                 grep ^PSC |
                 gzip -c > %s", vcf_file, "samp_stats.txt.gz"))


#### Read in intermediate files and Calculate Stats ####

## Read in the variants info and define het. frequency
vars <- fread("var_stats.txt.gz", header = FALSE, sep = "\t", na.strings = c("NA", "", "."))
names(vars) <- c("chrom", "pos", "dp", "maf", "f_miss", "ac_het", "ns")
vars$f_het <- vars$ac_het / vars$ns
vars$mean_nonmiss_dp <- vars$dp / vars$ns

## Read in the samples info and define both missing freq. and het. freq.
samps <- fread("samp_stats.txt.gz", header = FALSE, sep = "\t", 
               select = c(3, 4, 5, 6, 7, 8, 9, 10, 14))
names(samps) <- c("samp", "n_ref_hom", "n_nonref_hom", "n_het", "n_transt", 
                  "n_transv", "n_indel", "avg_depth", "n_miss")
samps$f_miss <- samps$n_miss/(samps$n_ref_hom + samps$n_nonref_hom + samps$n_het + samps$n_miss)
samps$f_het <- samps$n_het/(samps$n_ref_hom + samps$n_nonref_hom + samps$n_het)

## Calculate some variant-wise and sample-wise stats and write out
samp_stats <- describe(samps[, c("f_miss", "f_het", "avg_depth")], fast = TRUE)
vars_stats <- describe(vars[, c("dp", "mean_nonmiss_dp", "maf", "f_miss", "f_het")], fast = TRUE)
write.csv(samp_stats, "samplewise_statistics.csv")
write.csv(vars_stats, "variantwise_statistics.csv")


#### Plotting ####

## First plot sample-wise statistics
plot01 <- ggplot(samps, aes(x = f_miss)) +
  geom_histogram(binwidth = 0.01, color = "black", fill = "grey") +
  ggtitle("Sample-wise Missing Data Frequency") +
  xlab("Missing Data Freq.") +
  ylab("n Samples") +
  theme_minimal()
ggsave("plots/samplewise_missingness_hist.png", plot = plot01)

plot01 <- ggplot(samps, aes(x = f_het)) +
  geom_histogram(binwidth = 0.01, color = "black", fill = "grey") +
  ggtitle("Sample-wise Het. Call Frequency") +
  xlab("Het. Call Freq.") +
  ylab("n Samples") +
  theme_minimal()
ggsave("plots/samplewise_het_hist.png", plot = plot01)

## Now plot variant-wise statistics
plot01 <- ggplot(vars, aes(x = maf)) +
  geom_histogram(binwidth = 0.01, color = "black", fill = "grey") +
  ggtitle("Variant-wise MAF") +
  xlab("MAF") +
  ylab("n Variants") +
  theme_minimal()
ggsave("plots/MAF_hist.png", plot = plot01)

plot01 <- ggplot(vars, aes(x = f_miss)) +
  geom_histogram(binwidth = 0.01, color = "black", fill = "grey") +
  ggtitle("Variant-wise Missing Data Freq.") +
  xlab("Missing Data Freq.") +
  ylab("n Variants") +
  theme_minimal()
ggsave("plots/variantwise_missingness_hist.png", plot = plot01)

plot01 <- ggplot(vars, aes(x = f_het)) +
  geom_histogram(binwidth = 0.01, color = "black", fill = "grey") +
  ggtitle("Variant-wise Het. Call Freq.") +
  xlab("Het. Call Freq.") +
  ylab("n Variants") +
  theme_minimal()
ggsave("plots/variantwise_het_hist.png", plot = plot01)

plot01 <- ggplot(vars, aes(x = dp)) +
  geom_histogram(bins = 50, color = "black", fill = "grey") +
  ggtitle("Variant Total Depths") +
  xlab("Total Depth") +
  ylab("n Variants") +
  theme_minimal()
ggsave("plots/variantwise_total_depth_hist.png", plot = plot01)

plot01 <- ggplot(vars, aes(x = log10(dp))) +
  geom_histogram(bins = 50, color = "black", fill = "grey") +
  ggtitle("Log of Variant Total Depths") +
  xlab("Log10(Total Depth)") +
  ylab("n Variants") +
  theme_minimal()
ggsave("plots/variantwise_log_total_depth_hist.png", plot = plot01)

plot01 <- ggplot(vars, aes(x = mean_nonmiss_dp)) +
  geom_histogram(bins = 50, color = "black", fill = "grey") +
  ggtitle("Mean depth in samples with data") +
  xlab("Mean Depth") +
  ylab("n Variants") +
  theme_minimal()
ggsave("plots/variantwise_mean_nonmiss_depth_hist.png", plot = plot01)

## Remove intermediate files if specified 
if (delete_interfiles) {
  file.remove("samp_stats.txt.gz")
  file.remove("var_stats.txt.gz")
}
