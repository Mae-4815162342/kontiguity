#!/bin/bash
reference_genome=$1
result_path=$2
fastq_R1=$3
fastq_R2=$4
fastq=$5
hic_R1=$6
hic_R2=$7
enzymes=$8
hic_ref=$9
binning=$10
is_paired=$11
chromosome=$12

# 5. Hi-C mapping on new genome
# hicstuff pipeline call
hicstuff pipeline \
    --threads $proc \
    --enzyme $enzymes \
    --aligner bowtie2 \
    --no-cleanup \
    --filter \
    --duplicates \
    --prefix $hic_ref \
    --binning $binning \
    --matfmt bg2 \
    --plot \
    --force \
    --outdir $result_path/hic/$hic_ref \
    --genome $result_path/new_genome/new_genome \
    --distance-law \
    --mapping="iterative" \
    $hic_R1 \
    $hic_R2

# removing bam files as they stay in tmp with the --no-cleanup option
files=$(ls $result_path/hic/$hic_ref/tmp/*.bam)
for bam in $files;do
    rm "$bam"
done

# creating cool matrix
mkdir $result_path/hic/$hic_ref/images
cooler cload pairs \
    -c1 2 \
    -p1 3 \
    -c2 4 \
    -p2 5 \
    <(sed 1d $result_path/hic/$hic_ref/$hic_ref.chr.tsv | cut -f1-2):$binning \
    $result_path/hic/$hic_ref/tmp/$hic_ref.valid_idx_filtered.pairs $result_path/hic/$hic_ref/tmp/$hic_ref.valid_idx_filtered.pairs.cool

cooler zoomify \
    --nproc $proc \
    --resolutions 1000,2000,5000,10000,15000 \
    --balance \
    --out $result_path/hic/$hic_ref/tmp/$hic_ref.mcool \
    $result_path/hic/$hic_ref/tmp/$hic_ref.valid_idx_filtered.pairs.cool

for resolution in 1000 2000 5000 10000 15000
  do
    cooler show \
      --balanced \
      --dpi 500 \
      --out "$result_path/hic/$hic_ref/images/$ref.res$resolution.pdf" \
      --cmap="YlOrRd" \
      $result_path/hic/$hic_ref/tmp/$hic_ref.mcool::/resolutions/$resolution \
      $chromosome:0-
  done
