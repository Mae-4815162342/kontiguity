#!/bin/bash

outpath="test_results"
data_path="test_data/load"

kontiguity retrieve \
    -n S_cerevisiae \
    -o test_results \
    -i ${data_path}/S_cerevisiae \
    --wgs ${data_path}/wgs/FG0155_nxq_R1.fq.gz:${data_path}/wgs/FG0155_nxq_R2.fq.gz \
    --min-size 1000

# # external data
# kontiguity load \
#     -n "Vespula vulgaris" \
#     -o $outpath \
#     -i Vespula_vulgaris_1 \
#     --wgs ERR6054670,ERR6054671,ERR6054672

# # table test
# kontiguity load -n "table test" -o $outpath --table $data_path/test_dataset.csv

# # logan test
# kontiguity load \
#     -n "Vespula vulgaris" \
#     -o $outpath \
#     -i Vespula_vulgaris_1 \
#     --wgs ERR6054670,ERR6054671,ERR6054672 \
#     --logan