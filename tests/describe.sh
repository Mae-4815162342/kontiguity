#!/bin/bash

outpath="test_results"
data_path="test_data/describe"

## testing each statistic separately
# # Testing on known but un-described plasmid
# kontiguity describe \
#     -n E_histolytica \
#     -o test_results \
#     --mcool $data_path/Entamoeba_histolytica.mcool \
#     --description_binning 50000 \
#     --contigs AP023147.1 \
#     --index $data_path/E_histolytica \
#     --chromstart AP0 \
#     --min_chrom_size 1000000 \
#     --formats png,pdf

# table test
table_path=$data_path/parameter_table.csv
kontiguity describe \
    -o test_results \
    --chromstart AP0 \
    --contigs AP023147.1 \
    --table $table_path