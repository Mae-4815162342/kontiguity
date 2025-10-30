#!/bin/bash

outpath="test_results"
data_path="test_data/load"

# # local data
# bowtie2-build ${data_path}/S_cerevisiae.fna ${outpath}/S_cerevisiae
# kontiguity retrieve \
#     -n S_cerevisiae \
#     -o test_results \
#     -i ${data_path}/S_cerevisiae \
#     --wgs ${data_path}/wgs/FG0155_nxq_R1.fq.gz:${data_path}/wgs/FG0155_nxq_R2.fq.gz \
#     --min-size 1000 
    
# # external data
# kontiguity load \
#     -n "Vespula vulgaris" \
#     -o $outpath \
#     -r GCA_905475345.1 \
#     --wgs ERR6054670,ERR6054671,ERR6054672 \
#     --no_hic

# kontiguity retrieve \
#     -n "Vespula vulgaris" \
#     -o $outpath \
#     -r GCA_905475345.1 \
#     --wgs ERR6054670,ERR6054671,ERR6054672 \
#     --no_hic

# kontiguity retrieve \
#     -n "Vespula vulgaris" \
#     -o $outpath \
#     -i Vespula_vulgaris_1 \
#     --wgs ERR6054670,ERR6054671,ERR6054672 \
#     --no_tmp
#     --wgs ERR6054670,ERR6054671,ERR6054672 \
#     --no_tmp

# table test
tmp_out=/media/sardine/data_3/test_results
kontiguity retrieve -n "table test" -o $tmp_out --table $data_path/test_dataset.csv

# logan test
# logan test
# kontiguity load \
#     -n "Vespula vulgaris" \
#     -o $outpath \
#     -r GCA_905475345.1 \
#     --no_wgs \
#     --no_hic

# kontiguity retrieve \
#     -n "Vespula vulgaris" \
#     -o $outpath \
#     -i Vespula_vulgaris_1 \
#     --wgs ERR6054670,ERR6054671,ERR6054672 \
#     --logan \
#     --no_tmp