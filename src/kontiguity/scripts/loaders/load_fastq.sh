#!/bin/bash

fastq_ref=$1
outdir=$2
threads=$3

# TODO : put prefetch in a non-sbatch file as it will not work on the cluster
prefetch -O $outdir -X 100G $fastq_ref 2>$outdir/load_log.txt 1>$outdir/load_log.txt
fasterq-dump -O $outdir -e $threads -t $outdir $outdir/$fastq_ref 2>>$outdir/load_log.txt 1>>$outdir/load_log.txt