
# !/bin/bash

outpath=$1
accession=$2
process_stmp=$3 # id to inform the other processes the contigs have been loaded

local_path=$(realpath "$0")
local_dir="${local_path%/*}"
source "$local_dir/../lib/log.sh"

logs=$outpath/logs
if [ ! -d $outpath ];then
    mkdir $outpath
fi
if [ ! -d $logs ];then
    mkdir $logs
fi

log_info "[retrieve:$accession] Checking Logan for pre-built contigs"

## calls logan sevices on the required SRA accession number and retrieves contigs
aws s3 cp s3://logan-pub/c/${accession}/${accession}.contigs.fa.zst ${outpath}/contigs_${process_stmp}.fa.zst --no-sign-request 1>$logs/logan_log.txt 2>$logs/logan_log.txt
if [ -f ${outpath}/contigs_${process_stmp}.fa.zst ]; then
    zstd -d -f --rm -o ${outpath}/contigs_${process_stmp}.fa ${outpath}/contigs_${process_stmp}.fa.zst 1>>$logs/logan_log.txt 2>>$logs/logan_log.txt
    log_info "[retrieve:$accession] Contigs retrieved from Logan"
else
    log_info "[retrieve:$accession] Logan could not retrieve contigs, will build from scratch"
    echo Logan could not retrieve ${accession} contigs. Building from scratch. >> $logs/logan_log.txt
fi
