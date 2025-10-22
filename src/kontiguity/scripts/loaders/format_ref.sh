#!/bin/bash

fasta=$1
outfolder=$2
species=$3
sequence_types=$4

local_path=$(realpath "$0")
local_dir="${local_path%/*}"

# filtering chromosomes
if [ ! -f $outfolder/chromosomes.tsv ];then
    python3 $local_dir/format_ref.py $fasta $sequence_types $outfolder
else
    echo Fasta already formated for $species
fi

if [ ! -f $outfolder/$species.fa ];then
    # renaming fasta
    filtered_fasta=${outfolder}/$species.filtered.fa
    mv ${outfolder}/$species.filtered.fa ${outfolder}/$species.fa
fi

# creating bowtie index
bowtie2-build -q ${filtered_fasta%/*}/$species.fa ${filtered_fasta%/*}/$species