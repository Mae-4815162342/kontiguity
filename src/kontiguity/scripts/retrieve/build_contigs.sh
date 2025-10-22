#!/bin/bash

index=$1
outpath=$2
fastq_R1=$3
fastq_R2=$4
fastq=$5
is_paired=$6
min_len=$7
threads=$8
proc=$9

if [ -f ${outpath}/contigs.fa ];then
    echo Contigs already retrieved at ${outpath}/contigs.fa

else

    tmp_dir=$outpath/tmp
    if [ ! -d $outpath ];then
        mkdir $outpath
    fi
    if [ ! -d $tmp_dir ];then
        mkdir $tmp_dir
    fi

# 1. aligning assembly reads on reference genome and retrieving unmapped reads
    if [ "$is_paired" = true ]; then
        bowtie2 --un-conc $tmp_dir/unmapped.fastq -p $proc -x $index -1 $fastq_R1 -2 $fastq_R2 -S $tmp_dir/tmp_align.sam 1>>$outpath/contigs_log.txt 2>>$outpath/contigs_log.txt
    else
        bowtie2 --un $tmp_dir/unmapped.fastq -p $proc -x $index -U $fastq -S $tmp_dir/tmp_align.sam 1>>$outpath/contigs_log.txt 2>>$outpath/contigs_log.txt
    fi
    rm $result_path/tmp_align.sam

    # 2. assembly on unmapped reads
    # 2.1 delete one read if odd number of reads (for unpaired reads)
    if [ ! "$is_paired" = true ]; then
        lines=$(wc -l $tmp_dir/unmapped.fastq)
        lines=$(echo $lines | cut -d " " -f 1)
        nb_reads=$(($lines/4))
        if (( $nb_reads % 2 == 1))
        then 
            gawk -i inplace 'NR > 4' $tmp_dir/unmapped.fastq
        fi
    fi

    #2.2 genome assembly
    if [ "$is_paired" = true ]; then
        megahit --use-gpu -t $threads -1 $tmp_dir/unmapped.1.fastq -2 $tmp_dir/unmapped.2.fastq -o $tmp_dir --min-contig-len $min_len
    else
        megahit --use-gpu -t $threads --12 $tmp_dir/unmapped.fastq -o $tmp_dir --min-contig-len $min_len
    fi

    #2.3 retrieved contigs
    mv $tmp_dir/final.contigs.fa $outpath/contigs.fa

fi