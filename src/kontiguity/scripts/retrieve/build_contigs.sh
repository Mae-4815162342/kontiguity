#!/bin/bash

index=$1
outpath=$2
fastq_R1=$3
fastq_R2=$4
fastq=$5
is_paired=$6
min_len=$7
threads=$8
process_stmp=$9
no_tmp=${10}

logs=$outpath/logs
if [ ! -d $outpath ];then
    mkdir $outpath
fi
if [ ! -d $logs ];then
    mkdir $logs
fi

# writing contig retrieval informations
if [ "$is_paired" = "true" ]; then
    echo -e "Contigs retrieved from sequencing data:\n\t${fastq_R1}\n\t${fastq_R2}\nAligned on the reference genome indexed at ${index}\n" > $outpath/info.txt
else
    echo -e "Contigs retrieved from sequencing data:\n\t${fastq}\nAligned on the reference genome indexed at ${index}\n" \ > $outpath/info.txt
fi

if [ -f ${outpath}/contigs_${process_stmp}.fa ];then
    echo Contigs already retrieved at ${outpath}/contigs_${process_stmp}.fa > $logs/build_log.txt

    nb_sequences=$(cat $outpath/contigs_${process_stmp}.fa | grep -o ">" | wc -l)
    echo $nb_sequences contigs retrieved by Logan. >> $outpath/info.txt 
else

    tmp_dir=$outpath/tmp
    if [ ! -d $tmp_dir ];then
        mkdir $tmp_dir
    fi

    # 1. aligning assembly reads on reference genome and retrieving unmapped reads
    if [ "$is_paired" = "true" ]; then
        bowtie2 --un-conc $tmp_dir/unmapped.fastq -p $threads -x $index -1 $fastq_R1 -2 $fastq_R2 -S $tmp_dir/tmp_align.sam 1>$logs/build_log.txt 2>$logs/build_log.txt
        if [ "$no_tmp" = "true" ]; then
            rm $fastq_R1
            rm $fastq_R2
        fi
    else
        bowtie2 --un $tmp_dir/unmapped.fastq -p $threads -x $index -U $fastq -S $tmp_dir/tmp_align.sam 1>>$logs/build_log.txt 2>>$logs/build_log.txt
        if [ "$no_tmp" = "true" ]; then
            rm $fastq
        fi
    fi
    rm $tmp_dir/tmp_align.sam

    # 2. assembly on unmapped reads
    # 2.1 delete one read if odd number of reads (for unpaired reads)
    if [ ! "$is_paired" = "true" ]; then
        lines=$(wc -l $tmp_dir/unmapped.fastq)
        lines=$(echo $lines | cut -d " " -f 1)
        nb_reads=$(($lines/4))
        if (( $nb_reads % 2 == 1))
        then 
            gawk -i inplace 'NR > 4' $tmp_dir/unmapped.fastq
        fi
    fi

    #2.2 genome assembly
    tmp_dir_assembly=$outpath/assembly
    if [ -d $tmp_dir_assembly ];then
        rm -r $tmp_dir_assembly
    fi

    if [ "$is_paired" = "true" ]; then
        megahit -t $threads -1 $tmp_dir/unmapped.1.fastq -2 $tmp_dir/unmapped.2.fastq -o $tmp_dir_assembly --min-contig-len $min_len 1>>$logs/build_log.txt 2>>$logs/build_log.txt
    else
        megahit -t $threads --12 $tmp_dir/unmapped.fastq -o $tmp_dir_assembly --min-contig-len $min_len 1>>$logs/build_log.txt 2>>$logs/build_log.txt
    fi

    #2.3 retrieved contigs
    mv $tmp_dir_assembly/final.contigs.fa $outpath/contigs_${process_stmp}.fa

    if [ "$no_tmp" = "true" ]; then
        rm -r $tmp_dir_assembly
    fi

    nb_sequences=$(cat $outpath/contigs_${process_stmp}.fa | grep -o ">" | wc -l)
    echo $nb_sequences contigs retrieved by build. >> $outpath/info.txt 
fi