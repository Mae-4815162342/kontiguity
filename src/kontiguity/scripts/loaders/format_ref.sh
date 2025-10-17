#!/bin/bash

fasta=$1
outfolder=$2
species=$3
sequence_types=$4

local_path=$(realpath "$0")
local_dir="${local_path%/*}"

# filtering chromosomes
python3 $local_dir/format_ref.py $fasta $sequence_types $outfolder

# renaming fasta
filtered_fasta=${outfolder}/genome_filtered.fa
mv ${outfolder}/genome_filtered.fa ${outfolder}/$species.fa

# creating bowtie index
bowtie2-build ${filtered_fasta%/*}/$species.fa ${filtered_fasta%/*}/$species