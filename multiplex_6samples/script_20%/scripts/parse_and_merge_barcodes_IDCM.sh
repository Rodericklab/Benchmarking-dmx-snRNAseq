#!/bin/bash

# --------------------------------------------------
# Shell script to parse and merge cell barcode files
# --------------------------------------------------

# notes:
# - parse cell barcode files to contain unique sample IDs matching BAM files
# - there are 4 samples in IDCM dataset.

# runtime: seconds

# qsub -I -V -l nodes=1:ppn=3,pmem=60gb,walltime=01:00:00 -A lp_expcardiology -m bea -M yile.fu@kuleuven.be -d "$PWD"

# arguments:
# $1: directory for runtimes
# $2: directory for timestamp files
# $3: number of threads
# $4: output directory
# $5: outout directory for 6samples
# $6: sample ID 1
# $7: sample ID 2
# $8: sample ID 3
# $9: sample ID 4
# $10: sample ID 5
# $11: sample ID 6
# $12: short sample ID 1
# $13: short sample ID 2
# $14: short sample ID 3
# $15: short sample ID 4
# $16: short sample ID 5
# $17: short sample ID 6

# -----------------------------------
# start runtime
start=`date +%s`
# -----------------------------------

mkdir -p $5/barcodes_merged

# decompress cell barcode files for each sample
gunzip -c $4/$6/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > $5/barcodes_merged/barcodes_${12}.tsv
gunzip -c $4/$7/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > $5/barcodes_merged/barcodes_${13}.tsv
gunzip -c $4/$8/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > $5/barcodes_merged/barcodes_${14}.tsv
gunzip -c $4/$9/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > $5/barcodes_merged/barcodes_${15}.tsv
gunzip -c $4/${10}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > $5/barcodes_merged/barcodes_${16}.tsv
gunzip -c $4/${11}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz > $5/barcodes_merged/barcodes_${17}.tsv

# add unique sample IDs to cell barcodes for each sample
sed -i "s|\([A-Z]\+\)\-1|\1\-${12}|g" $5/barcodes_merged/barcodes_${12}.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-${13}|g" $5/barcodes_merged/barcodes_${13}.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-${14}|g" $5/barcodes_merged/barcodes_${14}.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-${15}|g" $5/barcodes_merged/barcodes_${15}.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-${16}|g" $5/barcodes_merged/barcodes_${16}.tsv
sed -i "s|\([A-Z]\+\)\-1|\1\-${17}|g" $5/barcodes_merged/barcodes_${17}.tsv

# merge files
cat $5/barcodes_merged/barcodes_${12}.tsv $5/barcodes_merged/barcodes_${13}.tsv \
    $5/barcodes_merged/barcodes_${14}.tsv $5/barcodes_merged/barcodes_${15}.tsv \
    $5/barcodes_merged/barcodes_${16}.tsv $5/barcodes_merged/barcodes_${17}.tsv > \
    $5/barcodes_merged/barcodes_merged.tsv


# -----------------------------------
# end runtime
end=`date +%s`
runtime=`expr $end - $start`

# save runtime
# mkdir -p $1/parse_and_merge_barcodes
echo runtime: $runtime seconds > $1/parse_and_merge_barcodes/runtime_parse_and_merge_barcodes_6samples.txt
# -----------------------------------


# -----------------------------------
# save timestamp file (for Snakemake)
# mkdir -p $2/parse_and_merge_barcodes
date > $2/parse_and_merge_barcodes/timestamp_parse_and_merge_barcodes_6samples.txt
# -----------------------------------

