#!/bin/bash

# -----------------------------------------------------------------------------------
# Shell wrapper for R script to generate awk lookup tables and updated barcodes files
# -----------------------------------------------------------------------------------

# Shell wrapper script to run R script "generate_awk_lookup_tables_doublets.R" 
# on computing cluster, required on our cluster since R module cannot be loaded 
# directly in Snakemake rule. If R is installed locally, this shell script is 
# not required, and the R script can be run directly instead.

# runtime: seconds

# qsub -I -V -l nodes=1:ppn=3,pmem=60gb,walltime=01:00:00 -A lp_expcardiology -m bea -M yile.fu@kuleuven.be -d "$PWD"

# arguments:
# $1: directory for runtimes
# $2: directory for timestamp files
# $3: output directory
# $4: scripts directory

Rscript $4/generate_awk_lookup_tables_doublets.R $1 $2 $3

