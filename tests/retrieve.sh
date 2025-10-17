#!/bin/bash

data_path=test_data/S_cerevisiae/dataset
kontiguity retrieve -n S_cerevisiae -o test_results -i ${data_path}/genomes/S_cerevisiae_1 --wgs ${data_path}/WGS/SRR35504920_1.fastq:${data_path}/WGS/SRR35504920_2.fastq --min-size 500