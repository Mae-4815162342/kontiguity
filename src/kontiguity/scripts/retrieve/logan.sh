
# !/bin/bash

outpath=$1
accession=$2
process_stmp=$3 # id to inform the other processes the contigs have been loaded


logs=$outpath/logs
if [ ! -d $outpath ];then
    mkdir $outpath
fi
if [ ! -d $logs ];then
    mkdir $logs
fi


## calls logan sevices on the required SRA accession number and retrieves contigs
aws s3 cp s3://logan-pub/c/${accession}/${accession}.contigs.fa.zst ${outpath}/contigs_${process_stmp}.fa.zst --no-sign-request 1>$logs/logan_log.txt 2>$logs/logan_log.txt
if [ -f ${outpath}/contigs_${process_stmp}.fa.zst ]; then
    zstd -d -f --rm -o ${outpath}/contigs_${process_stmp}.fa ${outpath}/contigs_${process_stmp}.fa.zst 1>>$logs/logan_log.txt 2>>$logs/logan_log.txt
else
    echo Logan could not retrieve ${accession} contigs. Building from scratch. >> $logs/logan_log.txt
fi
