#!/bin/bash

reference_genome=$1
result_path=$2
fastq_R1=$3
fastq_R2=$4
fastq=$5
is_paired=$6

proc=$(nproc --all) - 1

# 3. filtering assembly result by new alignment
bowtie2 -f --un $result_path/unmapped_assembly/unaligned_contigs.fa -p $proc -x $reference_genome -U  $result_path/unmapped_assembly/final.contigs.fa -S $result_path/unmapped_assembly/aligned_contigs.fa

# 4. building final reference genome
cat $reference_genome.fna $result_path/unmapped_assembly/unaligned_contigs.fa > $result_path/new_genome/new_genome.fna
bowtie2-build $result_path/new_genome/new_genome.fna $result_path/new_genome/new_genome