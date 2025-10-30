#!/bin/bash

reference_genome=$1
outpath=$2
min_len=$3
threads=$4
process_stmp=$5
no_tmp=$6

tmp_dir=$outpath/tmp
if [ ! -d $outpath ];then
    mkdir $outpath
fi
if [ ! -d $tmp_dir ];then
    mkdir $tmp_dir
fi

logs=$outpath/logs
if [ ! -d $logs ];then
    mkdir $logs
fi

local_path=$(realpath "$0")
local_dir="${local_path%/*}"

# filtering by sequence length
python3 $local_dir/filter_len.py $outpath/contigs_${process_stmp}.fa $min_len $tmp_dir/contigs_${process_stmp}.filtered.fa 1>$logs/filter_log.txt 2>$logs/filter_log.txt
mv $outpath/contigs_${process_stmp}.fa $tmp_dir/contigs_unfiltered.fa

nb_sequences=$(cat $tmp_dir/contigs_${process_stmp}.filtered.fa | grep -o ">" | wc -l)
echo $nb_sequences contigs kept with size \> $min_len kb. >> $outpath/info.txt

# filtering assembly result by new alignment
bowtie2 -f --un $tmp_dir/unaligned_contigs.fa -x $reference_genome -U $tmp_dir/contigs_${process_stmp}.filtered.fa -S $tmp_dir/aligned_contigs.fa 1>>$logs/filter_log.txt 2>>$logs/filter_log.txt

nb_sequences_unaligned=$(cat $tmp_dir/unaligned_contigs.fa | grep -o ">" | wc -l)
echo $nb_sequences_unaligned/$nb_sequences contigs kept after realignment. >> $outpath/info.txt

# building final reference genome
if [ -f "$reference_genome.fna" ]; then
    cp $reference_genome.fna $outpath/genome.fa 1>>$logs/filter_log.txt 2>>$logs/filter_log.txt
else
    cp $reference_genome.fa $outpath/genome.fa 1>>$logs/filter_log.txt 2>>$logs/filter_log.txt
fi
cat $tmp_dir/unaligned_contigs.fa >> $outpath/genome.fa

bowtie2-build -q $outpath/genome.fa $outpath/genome 1>>$logs/filter_log.txt 2>>$logs/filter_log.txt

mv $tmp_dir/unaligned_contigs.fa $outpath/contigs.fa 1>>$logs/filter_log.txt 2>>$logs/filter_log.txt

if [ "$no_tmp" = "true" ]; then
    rm -r $tmp_dir
fi
