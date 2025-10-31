#!/bin/bash

index=$1
outpath=$2
name=$3
hic_R1=$4
hic_R2=$5
enzymes=$6
binning=$7
threads=$8

tmp_dir=$outpath/tmp
if [ ! -d $outpath ];then
    mkdir $outpath
fi
if [ ! -d $tmp_dir ];then
    mkdir $tmp_dir
fi

logs=$outpath/logs
if [ ! -d $logs ];then
    mkdir $logs
fi

# calling hicstuff pipeline
hicstuff pipeline \
    --threads $threads \
    --force \
    --no-cleanup \
    --genome $index \
    --outdir $tmp_dir \
    --aligner bowtie2 \
    --enzyme $enzymes \
    --prefix $name \
    --filter \
    --duplicates \
    --matfmt bg2 \
    --distance-law \
    --mapping iterative \
    --binning $binning \
    --tmpdir $tmp_dir/hicstuff_tmp \
    $hic_R1 \
    $hic_R2 \
    1> $logs/mapping_log.txt 2> $logs/mapping_log.txt

echo -e "Map built from Hi-C data:\n\t${hic_R1}\n\t${hic_R1}\nAligned on the reference genome indexed at ${index}\n" > $outpath/info.txt

# cleaning up bams
files=$(ls $tmp_dir/hicstuff_tmp/*.bam)
for bam in $files;do
    rm "$bam"
done

# creating cool matrix
cooler cload pairs \
    -c1 2 \
    -p1 3 \
    -c2 4 \
    -p2 5 \
    <(sed 1d $tmp_dir/$name.chr.tsv | cut -f1-2):$binning \
    $tmp_dir/hicstuff_tmp/$name.valid_idx_filtered.pairs $outpath/$name.$binning.cool \
    1>> $logs/mapping_log.txt 2>> $logs/mapping_log.txt

echo Initial cool matrix built for binning $binning bp at $outpath/$name.$binning.cool. >> $outpath/info.txt

# balancing
cooler balance \
    --nproc $threads \
    --force \
    $outpath/$name.$binning.cool \
    1>> $logs/mapping_log.txt 2>> $logs/mapping_log.txt

echo $outpath/$name.$binning.cool has been balanced. >> $outpath/info.txt

mv $tmp_dir/$name.chr.tsv $outpath/$name.chr.tsv 
mv $tmp_dir/$name.frags.tsv $outpath/$name.frags.tsv 
mv $tmp_dir/$name.distance_law.txt $outpath/$name.distance_law.txt
