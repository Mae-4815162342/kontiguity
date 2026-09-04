#!/bin/bash

outpath=$1
name=$2
original_binning=$3
binnings=$4
threads=$5
no_tmp=$6

local_path=$(realpath "$0")
local_dir="${local_path%/*}"
source "$local_dir/../lib/log.sh"

if [ ! -d $outpath ];then
    mkdir $outpath
fi

tmp_dir=$outpath/tmp

logs=$outpath/logs
if [ ! -d $logs ];then
    mkdir $logs
fi

for binning in $binnings; do

    log_info "[map:$name] Rebinning $name.$original_binning.cool at ${binning}bp"
    echo Rebinning $name.$original_binning.cool at $binning bp. >> $logs/mapping_log.txt

    cooler cload pairs \
        -c1 2 \
        -p1 3 \
        -c2 4 \
        -p2 5 \
        <(sed 1d $outpath/$name.chr.tsv | cut -f1-2):$binning \
        $tmp_dir/hicstuff_tmp/$name.valid_idx_filtered.pairs $outpath/$name.$binning.cool \
        1>> $logs/mapping_log.txt 2>> $logs/mapping_log.txt

    # balancing
    cooler balance \
        --nproc $threads \
        --force \
        $outpath/$name.$binning.cool \
        1>> $logs/mapping_log.txt 2>> $logs/mapping_log.txt

done

echo Cool matrix built and balanced for each binning "(${binnings})". >> $outpath/info.txt
log_info "[map:$name] Rebinning complete for ($binnings)"

if [ "$no_tmp" = "true" ]; then
    rm -r $tmp_dir
fi
