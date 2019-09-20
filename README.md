# hpc_seq_process

Scripts for processing sequencing data on HPC clusters

## Description

This repository is a collection of scripts for processing sequencing data
(currently just DNA) on high-performance computing (HPC) clusters. It is
currently geared towards clusters using SLURM job schedulers.

## Parallel execution

Many scripts in this repository end with "\_parallel.sh". These scripts are
intended to work on files corresponding to a single sample. They may be used
with a dispatcher script, which coordinates parallel or sequential execution
for multiple samples. These scripts are located in code/parallel_dispatch.
Generally, on clusters arrayer.sh is the preferred script for this purpose.

All of the dispatching scripts simply supply integers within a range starting at
one. This range of integers is then used to subset lines out of a text file,
listing each sample (the path to this file is defined in whatever script is
being controlled).

## Execution on local machines

The code in this repository can be run outside of a HPC setting. In this
use case, parallelizer.sh or looper.sh would be used for parallel (or sequential)
dispatch. Lines within scripts which start with "module" must be commented out.
These lines also show the dependencies for each script. Outside of an HPC
setting, it is usually easiest to use Bioconda environments 
(https://bioconda.github.io/user/install.html) to handle dependencies.

## Typical Workflow - Short Variant Identification

Steps below would typically be run using a dispatcher script for parallel execution.

For short variant identification, a common pre-processing workflow would be:

1. run code/concat_fastqs.sh to generate a single mated pair of fastq files for
each individual sample
2. run code/fastq_filt_trim/bbduk_filt_trim_[PE/SE]_parallel.sh to clean fastq reads
for paired-end or single-end reads, respectively
3. run code/alignment/bowtie2_[PE/SE]_align_parallel.sh to align reads for
paired end or single end reads using bowtie2 and create single sample bam alignment files
4. run code/filt_bam_files_parallel.sh to filter bam files after alignment (this
step is optional, as bcftools mpileup can perform filtering during variant calling)

For variant calling, the following steps may then be performed after bam file
generation:

1. run code/call_variants/call_variants.sh
2. run code/call_variants/filter_raw_vcf.sh
3. run code/call_variants/rename_annotate_split_vcf.sh to (optionally) rename samples,
predict variant transcriptional effects using a genome annotation file, and create subset 
VCF files consisting of just indels and just SNPs.
