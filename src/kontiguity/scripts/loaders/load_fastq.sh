#!/bin/bash

fastq_ref=$1
outdir=$2
threads=$3

fasterq-dump -O $outdir -e $threads -t $outdir/tmp $fastq_ref