#!/bin/bash

fastq_ref=$1
outdir=$2
threads=$3

local_path=$(realpath "$0")
local_dir="${local_path%/*}"
source "$local_dir/../lib/log.sh"

mkdir -p "$outdir"

log_info "Loading fastq(s) for accession $fastq_ref"

# Ask ENA which fastq.gz files exist for this accession (paired-end runs
# return two, single-end runs return one, separated by ';').
urls=$(wget -q -O - \
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${fastq_ref}&result=read_run&fields=fastq_ftp&format=tsv" \
    | tail -n +2 | cut -f2) 2>>$outdir/load_log.txt 1>>$outdir/load_log.txt

if [ -z "$urls" ]; then
    log_error "${fastq_ref}: no fastq_ftp entries returned by ENA (accession may not be mirrored there yet)"
    echo "!!! ${fastq_ref}: no fastq_ftp entries returned by ENA (accession may not be mirrored there yet)" 2>>$outdir/load_log.txt 1>>$outdir/load_log.txt
    exit 1
fi

echo $urls 2>>$outdir/load_log.txt 1>>$outdir/load_log.txt
IFS=';' read -ra url_array <<< "$urls"
for url in "${url_array[@]}"; do
    fname=$(basename "$url")
    log_info "Downloading ${fname} for ${fastq_ref}"
    echo ">>> Downloading ${fname}" 2>>$outdir/load_log.txt 1>>$outdir/load_log.txt
    wget -c -q --tries=5 --timeout=60 -O "${outdir}/${fname}" "https://${url}" 2>>$outdir/load_log.txt 1>>$outdir/load_log.txt
done

log_info "Done loading fastq(s) for ${fastq_ref}, files are in ${outdir}"
echo "Done. FASTQ files are in ${outdir}" 2>>$outdir/load_log.txt 1>>$outdir/load_log.txt
