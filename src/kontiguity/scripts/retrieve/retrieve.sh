#!/bin/bash

reference_genome=$1
result_path=$2
fastq_R1=$3
fastq_R2=$4
fastq=$5
is_paired=$6

proc=$(nproc --all) - 1

# 1. aligning assembly reads on reference genome and retrieving unmapped reads
if [ ! -d $result_path ];then
    mkdir $result_path
fi

if [ "$is_paired" = true ]; then
    bowtie2 --un-conc $result_path/unmapped.fastq -p $proc -x $reference_genome -1 $fastq_R1 -2 $fastq_R2 -S $result_path/tmp_align.sam
else
    bowtie2 --un $result_path/unmapped.fastq -p $proc -x $reference_genome -U $fastq -S $result_path/tmp_align.sam
fi
rm $result_path/tmp_align.sam

# 2. assembly on unmapped reads

# 2.1 delete one read if odd number of reads (for unpaired reads)
if [ ! "$is_paired" = true ]; then
    lines=$(wc -l $result_path/unmapped.fastq)
    lines=$(echo $lines | cut -d " " -f 1)
    nb_reads=$(($lines/4))
    if (( $nb_reads % 2 == 1))
    then 
        gawk -i inplace 'NR > 4' $result_path/unmapped.fastq
    fi
fi

 #2.2 genome assembly
if [ -d "${result_path}/unmapped_assembly" ]; then
    rm -rf $result_path/unmapped_assembly
fi
if [ "$is_paired" = true ]; then
    megahit --use-gpu -t $proc -1 $result_path/unmapped.1.fastq -2 $result_path/unmapped.2.fastq -o $result_path/unmapped_assembly --min-contig-len $binning
else
    megahit --use-gpu -t $proc --12 $result_path/unmapped.fastq -o $result_path/unmapped_assembly --min-contig-len $binning
fi

if [ ! -f "${result_path}/new_genome" ]; then
    mkdir $result_path/new_genome
fi

# 3. filtering assembly result by new alignment
bowtie2 -f --un $result_path/unmapped_assembly/unaligned_contigs.fa -p $proc -x $reference_genome -U  $result_path/unmapped_assembly/final.contigs.fa -S $result_path/unmapped_assembly/aligned_contigs.fa

# 4. building final reference genome
cat $reference_genome.fna $result_path/unmapped_assembly/unaligned_contigs.fa > $result_path/new_genome/new_genome.fna
bowtie2-build $result_path/new_genome/new_genome.fna $result_path/new_genome/new_genome