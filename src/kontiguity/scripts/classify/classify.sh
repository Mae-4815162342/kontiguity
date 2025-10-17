#!/bin/bash

### passing matrix, new contigs fasta and original genome fasta file and output path
python3 $installation_path/pipeline_classification/classify.py \
    $result_path/hic/$hic_ref/tmp/$hic_ref.mcool \
    $reference_genome.fna \
    $result_path/classification