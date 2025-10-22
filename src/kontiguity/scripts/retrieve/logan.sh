
# !/bin/bash

outpath=$1
accession=$2

## calls logan sevices on the required SRA accession number and retrieves contigs
wget -nv -O ${outpath}/contigs.fa.zst https://s3.amazoneaws.com/logan-pub/c/${accession}/${accession}.contigs.fa.zst
zstd -d -o ${outpath}/contigs.fa ${outpath}/contigs.fa.zst