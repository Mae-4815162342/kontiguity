#!/bin/bash

outpath="test_results"
data_path="test_data/load"

# local data
# bowtie2-build ${data_path}/S_cerevisiae.fna ${data_path}/S_cerevisiae
kontiguity map \
    -n S_cerevisiae \
    -o test_results \
    -i ${data_path}/S_cerevisiae \
    --hic $data_path/hic/FG0153_nxq_R1.fq.gz:$data_path/hic/FG0153_nxq_R2.fq.gz \
    --enzymes HinfI,DpnII \
    --format mcool \
    --binnings 10000,15000 \
    --no_tmp

# # external data
# kontiguity load \
#     -n "Vespula vulgaris" \
#     -o $outpath \
#     -r GCA_905475345.1 \
#     --hic ERR6054673,ERR6054674,ERR6054675 \
#     --no_wgs

# kontiguity map \
#     -n "Vespula vulgaris" \
#     -o $outpath \
#     -i Vespula_vulgaris_1 \
#     --hic ERR6054673,ERR6054674,ERR6054675 \
#     --no_tmp

# # table test
# tmp_out=/media/sardine/data_3/test_results
# kontiguity map -n "table test" -o $tmp_out --table $data_path/test_dataset.csv