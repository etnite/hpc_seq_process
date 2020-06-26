# hpc_seq_process

Scripts for processing sequencing data on HPC clusters

Author: [Brian Ward](https://brianpward.net/)  
Email: [brian@brianpward.net](mailto:brian@brianpward.net)  
Github: [etnite](https://github.com/etnite)  
License: [GPLv3](https://opensource.org/licenses/GPL-3.0)

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

This repository is primarily set up to use the samtools/bcftools variant calling pipeline. 
Steps below would typically be run using a dispatcher script for parallel execution. 
Note that for samtools/bcftools variant calling, read groups (e.g. different
flowcell/lane combinations) can be concatenated into single fastq files for each sample.
However, for modern GATK workflows, it is more common to work with separate files for
individual read groups.

For short variant identification, a common pre-processing workflow would be:

1. Run code/fastq_filt_trim/bbduk_filt_trim_[PE/SE]_parallel.sh to clean fastq reads
for paired-end or single-end reads, respectively
2. Run code/alignment/bowtie2_[PE/SE]_align_parallel.sh to align reads for
paired end or single end reads using bowtie2 and create single sample bam alignment files
3. (optional) run code/bam_manipulations/merge_bams_parallel.sh to merge and filter BAM files (e.g.
to create one filtered BAM file per sample)
5. (optional) - for some applications, e.g. practical haplotype graph creation, convert BAM files
into single-sample gVCFs using code/gvcf_creation/create_single_samp_gvcf_parallel.sh

For variant calling, the following steps may then be performed after bam file
generation:

1. Run code/call_variants/mpileup_call_parallel.sh - this script can run in parallel on multiple
regions of the genome, for instance those created with code/make_regions_file.py
2. Run code/call_variants/concat_sort_vcfs.sh to concatenate individual region VCFs/BCFs into one file
3. (optional) Run code/call_variants/rename_annotate_split_vcf.sh to rename samples,
predict variant transcriptional effects using a genome annotation file, and create subset 
VCF files consisting of just indels and just SNPs.
