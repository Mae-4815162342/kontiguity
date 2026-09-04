#!/bin/bash

fasta=$1
outfolder=$2
species=$3
sequence_types=$4

local_path=$(realpath "$0")
local_dir="${local_path%/*}"
source "$local_dir/../lib/log.sh"

# filtering chromosomes
if [ ! -f $outfolder/chromosomes.tsv ];then
    log_info "Formatting reference for $species"
    python3 $local_dir/format_ref.py $fasta $species $sequence_types $outfolder
else
    log_info "Fasta already formated for $species"
fi

if [ ! -f $outfolder/$species.fa ];then
    # renaming fasta
    filtered_fasta=${outfolder}/$species.filtered.fa
    mv ${outfolder}/$species.filtered.fa ${outfolder}/$species.fa
fi

# creating bowtie index
log_info "Building bowtie2 index for $species"
bowtie2-build -q ${outfolder}/$species.fa ${outfolder}/$species
