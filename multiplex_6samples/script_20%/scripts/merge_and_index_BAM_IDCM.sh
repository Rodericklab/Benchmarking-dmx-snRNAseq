#!/bin/bash

# ------------------------------------------------
# Shell script to merge and index parsed BAM files
# ------------------------------------------------

# notes:
# - BAM files from Cell Ranger are already position sorted, so do not need to sort
# - there are 4 samples in IDCM dataset.

# runtime: ~4-6 hours

# qsub -I -V -l nodes=1:ppn=3,pmem=60gb,walltime=06:00:00 -A lp_expcardiology -m bea -M yile.fu@kuleuven.be -d "$PWD"

# arguments:
# $1: directory for runtimes
# $2: directory for timestamp files
# $3: number of threads
# $4: output directory
# $5: output directory for 6 samples
# $6: sample ID 1
# $7: sample ID 2
# $8: sample ID 3
# $9: sample ID 4
# $10: sample ID 5
# $11: sample ID 6
# -----------------------------------
# start runtime
start=`date +%s`
# -----------------------------------

mkdir -p $5/bam_merged

# merge BAM files
samtools merge $5/bam_merged/bam_merged.bam \
$4/$6/outs/possorted_genome_bam_parsed.bam $4/$7/outs/possorted_genome_bam_parsed.bam $4/$8/outs/possorted_genome_bam_parsed.bam $4/$9/outs/possorted_genome_bam_parsed.bam $4/${10}/outs/possorted_genome_bam_parsed.bam $4/${11}/outs/possorted_genome_bam_parsed.bam

# index merged BAM
samtools index $5/bam_merged/bam_merged.bam


# -----------------------------------
# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
# mkdir -p $1/merge_and_index_BAM
echo runtime: $runtime seconds > $1/merge_and_index_BAM/runtime_merge_and_index_BAM_6samples.txt
# -----------------------------------


# -----------------------------------
# save timestamp file (for Snakemake)
# mkdir -p $2/merge_and_index_BAM
date > $2/merge_and_index_BAM/timestamp_merge_and_index_BAM_6samples.txt
# -----------------------------------

