#!/bin/bash

outpath="test_results"
data_path="test_data/describe"

## testing each statistic separately
# Testing on known but un-described plasmid
kontiguity describe \
    -n E_histolytica \
    -o test_results \
    --mcool $data_path/Entamoeba_histolytica.mcool \
    --binning 50000 \
    --contigs AP023147.1 \
    --fasta $data_path/E_histolytica.fna \
    --chromstart AP0 \
    --min_chrom_size 1000000 \
    --formats png,pdf