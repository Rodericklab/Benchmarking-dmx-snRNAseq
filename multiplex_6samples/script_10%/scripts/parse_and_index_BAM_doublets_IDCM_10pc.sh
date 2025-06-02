#!/bin/bash

#########################################
# Shell script to run doublets simulation
#########################################

# This script runs part of the doublets simulation by:
# (i) converting BAM to SAM
# (ii) parsing through the merged SAM file to replace and combine some 
# percentage of cell barcodes
# (iii) converting SAM back to BAM

# The modifed BAM file can then be used by cellSNP and Vireo (or alternative 
# tools) in the following scripts.

# Notes:
# - lookup table to use in awk command (and updated list of cell barcodes) are 
# generated in previous step with script "generate_awk_lookup_tables_doublets.R"
# - for more details on how to use awk for multiple replacements see: 
# https://stackoverflow.com/questions/14234907/replacing-values-in-large-table-using-conversion-table

# runtime: ~1 day

# qsub -I -V -l nodes=1:ppn=3,pmem=60gb,walltime=12:00:00 -A lp_expcardiology -m bea -M yile.fu@kuleuven.be -d "$PWD"

# arguments:
# $1: directory for runtimes
# $2: directory for timestamp files
# $3: number of threads
# $4: output directory
# $5: output directory for 6 samples
# $6: percentage of doublets for simulation scenario (formatted as e.g. "10pc")


# -----------------------------------
# start runtime
start=`date +%s`
# -----------------------------------


# -----------------------------
# Parse through merged BAM file
# -----------------------------

# parse through merged BAM file to combine some percentage of cell barcodes
# for one simulation scenario (dataset, percentage of doublets)

# note hyphen for argument order
samtools view -h $5/bam_merged/bam_merged.bam | \
awk \
'NR==1 { next } FNR==NR { a[$1]=$2; next } (i=gensub(/.*CB\:Z\:([A-Za-z]+\-[A-Za-z0-9]+).*/, "\\1", 1, $0)) in a { gsub(i, a[i]) }1' \
$5/doublets_sims/$6/lookup_table_doublets_IDCM_$6.tsv - | \
samtools view -bo $5/doublets_sims/$6/bam_merged_doublets_IDCM_$6.bam


# ---------
# Index BAM
# ---------

samtools index $5/doublets_sims/$6/bam_merged_doublets_IDCM_$6.bam


# -----------------------------------
# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
# mkdir -p $1/$5/doublets_sims/$6/parse_and_index_BAM_doublets
echo runtime: $runtime seconds > $1/IDCM/doublets_sims/$6/parse_and_index_BAM_doublets/runtime_parse_and_index_BAM_doublets_IDCM_{$6}_6samples.txt
# -----------------------------------


# -----------------------------------
# save timestamp file (for Snakemake)
# mkdir -p $2/$5/doublets_sims/$6/parse_and_index_BAM_doublets
date > $2/IDCM/doublets_sims/$6/parse_and_index_BAM_doublets/timestamp_parse_and_index_BAM_doublets_IDCM_{$6}_6samples.txt
# -----------------------------------

