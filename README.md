# Kontiguity: tool for eukaryotes contigs retrieval from genomic data and classification

Kontiguity is a python and bash pipeline created to retrieve unidentified contigs from eukaryotes genomic data, and to classify said contigs based on genomic contact data (Hi-C).

## Installation

```bash
pip install kontiguity
```

For development:

```bash
git clone https://github.com/Mae-4815162342/kontiguity.git
cd kontiguity
pip install -e .
```

## Presenting pipeline

Kontiguity is based on a pipeline of four subfunctions:

- **load** which serves data retrieval and formating, and provides a scrapping method to build a dataset from DToL (ref to add).

- **retrieve** for the retrieval of new contigs from WGS reads aligned on a reference genome.

- **contact-map** maps Hi-C data on new genomes (using hicstuff), building mcool files.

- **classify** for the contigs contacts classification, based on a plasmid-detection-oriented model at this day (a larger model can be provided later).

Those four functions can be called individually or in order with the **pipeline** command, which provides an option to start at any step.

## Usage

### Loading a dataset

```bash
kontiguity load -n Saccaromyces_cerevisiae -o outfolder -r S_cerevisiae.fa --chroms chromosome.tsv --wgs fastq_wgs.fq.gz --hic fastq_R1.fq.gz,fastq_R2.fq.gz

kontiguity load -o outfolder --table samples.csv

kontiguity load -o outfolder --dtol
```

Options:
```
-n/--name       name of the experiment (recommanded: species name. info: spaces are not allowed and will be replaced by _.)
-o/--outpath    output folder path, created if non-existent
-r/--ref        path to the reference genome fasta OR the GCA reference which will automatically be loaded from ENA database.
--chroms        path to a chromosome information file detailing the type of each sequence present in the reference (Mandatory column heads: ["id", "sequence_type", "sequence_name"]). "sequence_type" must be in the ENA database format : ["chromosome", "organelle", ...]. Required only for a local fasta, GCA referenced genomes will have the chromosome.tsv generated.
--wgs           path to the WGS fastq(s) OR SRA accession. If paired and local, provide both fastqs comma-separated.
--hic           path to the Hi-C fastq(s) OR SRA accession. If paired and local, provide both fastqs comma-separated.
--table         path to a csv table providing the data parameters (Mandatory column heads: ["name", "ref", "wgs", "hic"]). 
--dtol          if selected, a data table will be created and loaded from the Darwin Tree of Life project database.
```

Output:
At the outfolder/name location, the following file arborescence is built:
```
outfolder/name
    └── dataset
        ├── genomes
        │   ├── genome1
        |   |   ├── *bowtie index*
        |   |   └── chromosome.tsv *chromosomes informations*
        |   └── ...
        └── fastqs
            ├── WGS
            │   ├── *downloaded fastq files*
            |   └── summup.csv *summup of the paths to the fastq files*
            └── HiC
                ├── *downloaded fastq files*
                └── summup.csv *summup of the paths to the fastq files*
```
If several species are provided, each gets its individual arboresence.

### Retrieving contigs

This command will retrieve the contigs from a WGS aligned on the reference genome by assembling the unaligned reads. The new contigs are added at the end of the reference genome and a new bowtie index is generated.

```bash
kontiguity retrieve -n Saccaromyces_cerevisiae -o outfolder -i S_cerevisiae --wgs fastq_wgs.fq.gz

kontiguity retrieve -n Saccaromyces_cerevisiae -o outfolder --min-size 1500 --table summup.csv
```

Options:
```
-n/--name       name of the experiment (recommanded: species name. info: spaces are not allowed and will be replaced by _.)
-o/--outpath    output folder path, created if non-existent.
-i/--index      path to the reference genome index.
--min-size      minimum size of the kept contigs in bp (dflt: 1000).
--wgs           path to the WGS fastq(s).
--table         path to a csv table providing the data parameters (Mandatory column heads: ["name", "index", "wgs"]). 
```

Output:
At the outfolder/name location, the following file arborescence is built:
```
outfolder/name
    └── contigs
```

### Mapping Hi-C

Call to hicue to map each provided Hi-C fastqs on each provided genome. In the pipeline, the maps are generated on the new genomes with retrieved contigs from the **retrieve** command.

```bash
kontiguity map -n Saccaromyces_cerevisiae -o outfolder -g genome_path --hic fastq_R1.fq.gz,fastq_R2.fq.gz --enzymes DpnII,HinfI
```

Options:
```
-n/--name       name of the experiment (recommanded: species name. info: spaces are not allowed and will be replaced by _).
-o/--outpath    output folder path, created if non-existent.
-i/--index      path to the genome index (from retrieved in the pipeline)
--hic           path to the Hi-C fastq(s).
--enzymes       Hi-C restriction enzymes (dflt: DpnII,HinfI). The default enzymes where chosen in regard of the Arima Hi-C kit (ref).
--table         path to a csv table providing the data parameters (Mandatory column heads: ["name", "index", "hic", "enzymes"]). 
+ hicstuff parameters (ref)
```

Output:
At the outfolder/name location, the following file arborescence is built:
```
outfolder/name
    └── hic
```

### Classifying contigs contacts

```bash
kontiguity classify -n Saccaromyces_cerevisiae -o outfolder --chroms chromosome.tsv --mcool S_cerevisiae_on_genome1.mcool --model regression.hdf5 
```

Options:
```bash
```

Output:
At the outfolder/name location, the following file arborescence is built:
```
outfolder/name
    └── classification
```

### Pipeline

```bash
kontiguity pipeline -n Saccaromyces_cerevisiae -o outfolder -r S_cerevisiae.fa --wgs fastq_wgs.fq.gz --hic fastq_R1.fq.gz,fastq_R2.fq.gz
```

Options:
```bash
```

TODO re-write with real test data

## Outputs

TODO

## Classification model

TODO