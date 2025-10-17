#!/bin/bash

outpath="test_results"
data_path="test_data/load"
echo $data_path
# # local data
# kontiguity load \
#     -n "Saccharomyces cerevisiae" \
#     -o $outpath \
#     -r $data_path/S_cerevisiae.fna \
#     --chroms $data_path/chromosomes.csv \
#     --wgs $data_path/wgs/FG0155_nxq_R1.fq.gz:$data_path/wgs/FG0155_nxq_R2.fq.gz \
#     --hic $data_path/hic/FG0153_nxq_R1.fq.gz:$data_path/hic/FG0153_nxq_R2.fq.gz


# external data
kontiguity load \
    -n "Vespula vulgaris" \
    -o $outpath \
    -r GCA_905475345.1 \
    --wgs ERR6054670,ERR6054671,ERR6054672 \
    --hic ERR6054673,ERR6054674,ERR6054675

# # table test
# kontiguity load -n "table test" -o $outpath --table $data_path/test_dataset.csv

# # dtol test
# kontiguity load -n dtol_test -o $outpath --dtol