#!/bin/bash

fasta=$1
output_folder=$2
species=$3

local_path=$(realpath "$0")
local_dir="${local_path%/*}"
source "$local_dir/../lib/log.sh"

if [ ! -d $output_folder ]; then
    mkdir $output_folder
fi

if [ -f $output_folder/$species.all_seqs.fa ];then
    log_info "Fasta already retrieved for $species"
else
    log_info "Retrieving fasta for $species from $fasta"
    HTTP_CODE=$(curl --silent --output $output_folder/$species.all_seqs.fa.gz --write-out "%{http_code}" $fasta)
    if [[ ${HTTP_CODE} -lt 200 || ${HTTP_CODE} -gt 299 ]] ; then
        log_error "Genome not found for ${species}: error $HTTP_CODE"
    else
        gunzip -f $output_folder/$species.all_seqs.fa.gz
        log_info "Fasta retrieved for $species"
    fi
fi
