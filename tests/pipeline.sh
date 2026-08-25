#!/bin/bash

outpath="test_results/pipeline"
data_path="test_data/load"

# Testing the pipeline
kontiguity pipeline \
    -n "Saccharomyces cerevisiae" \
    -o $outpath \
    -i $data_path/S_cerevisiae \
    --chroms $data_path/chromosomes.csv \
    --wgs $data_path/wgs/FG0155_nxq_R1.fq.gz:$data_path/wgs/FG0155_nxq_R2.fq.gz \
    --hic $data_path/hic/FG0153_nxq_R1.fq.gz:$data_path/hic/FG0153_nxq_R2.fq.gz \
    --first_step describe \
    --formats png,pdf \
    --binnings 10000,50000
