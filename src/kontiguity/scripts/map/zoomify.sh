#!/bin/bash

outpath=$1
name=$2
binning=$3
zoomings=$4
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

log_info "[map:$name] Building mcool with resolutions $zoomings"

# creating mcool
cooler zoomify \
    --nproc $threads \
    --resolutions $zoomings \
    --balance \
    --out $outpath/$name.mcool \
    $outpath/$name.$binning.cool \
    1>> $logs/mapping_log.txt 2>> $logs/mapping_log.txt

mv $outpath/$name.$binning.cool $tmp_dir/$name.$binning.cool

echo Mcool matrix built\; contains cool matrix built and balanced for each binning "(${zoomings})". >> $outpath/info.txt
log_info "[map:$name] mcool built at $outpath/$name.mcool"

if [ "$no_tmp" = "true" ]; then
    rm -r $tmp_dir
fi
