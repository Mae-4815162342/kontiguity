#!/bin/bash

fasta=$1
output_folder=$2
species=$3

if [ ! -d $output_folder ]; then
    mkdir $output_folder
fi

HTTP_CODE=$(curl --silent --output $output_folder/genome.fa.gz --write-out "%{http_code}" $fasta)
if [[ ${HTTP_CODE} -lt 200 || ${HTTP_CODE} -gt 299 ]] ; then
    echo Genome not found for ${species}: error $HTTP_CODE
else
    gunzip -f $output_folder/genome.fa.gz
    echo Fasta retrieved for $species
fi
