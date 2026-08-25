import click
from .types import *
import kontiguity.pipeline as kpipeline

@click.command("pipeline")
@click.option(
    '-n',
    '--name', 
    type=str,
    default="",
    help="name of the experiment (recommanded: species name. info: spaces are not allowed and will be replaced by _.)"
)
@click.option(
    '-o',
    '--outpath',
    type=str,
    help="output folder path, created if non-existent"
)
@click.option(
    "--first_step",
    default = "load",
    type=click.Choice(["load", "retrieve", "map", "describe"]),
    help='starting point of the pipeline. Requires the arboresence of any previous step to exist (dflt: load).'
)
@click.option(
    "--last_step",
    default = "describe",
    type=click.Choice(["load", "retrieve", "map", "describe"]),
    help='stopping point of the pipeline (dflt: describe).'
)
### data parameters
@click.option(
    "--wgs",
    type=PAIR_LIST,
    help="comma-separated list of paths to the WGS fastq(s) OR SRA accession. If paired and local, provide both fastqs separated by : ."
)
@click.option(
    "--hic",
    type=PAIR_LIST,
    help="comma-separated list of paths to the Hi-C fastq(s) OR SRA accession. If paired and local, provide both fastqs separated by : "
)
@click.option(
    "-i",
    "--index",
    type=str,
    help="path to the reference genome index."
)
@click.option(
    "--chroms",
    type=str,
    help='path to a chromosome information file detailing the type of each sequence present in the reference (Mandatory column heads: ["id", "sequence_type", "sequence_name"]). "sequence_type" must be in the ENA database format : ["chromosome", "organelle", ...]. Required only for a local fasta, GCA referenced genomes will have the chromosome.tsv generated.'
)
@click.option(
    "--table",
    type=str,
    help='path to a csv table providing the data parameters: \n for load: (Mandatory column heads: ["name", "ref", "wgs", "hic"]). See test_data/load/test_dataset.csv for format example.'
)
### load options
@click.option(
    "-r",
    "--ref",
    type=str,
    help="path to the reference genome fasta OR the GCA reference which will automatically be loaded from ENA database."
)
@click.option(
    "--dtol",
    is_flag=True,
    default=False,
    help="if selected, a data table will be created and loaded from the Darwin Tree of Life project database."
)
@click.option(
    "--no_wgs",
    is_flag=True,
    default=False,
    help="if selected, the WGS data is not loaded nor required. Usefull for calling Logan in contigs reconstruction (see retrieve command)."
)
@click.option(
    "--no_hic",
    is_flag=True,
    default=False,
    help="if selected, the Hi-C data is not loaded nor required. Usefull when using kontiguity for contigs retrieval only."
)
### retrieve options
@click.option(
    "--min-size",
    type=int,
    default=1000,
    help='minimum size of the kept contigs in bp.'
)
@click.option(
    "--logan",
    is_flag=True,
    default=False,
    help="if selected, will call to the AWS Logan database from SRA accession number to retrieve contigs. If contigs are found, they will not be built from scratch. (dflt: False)"
)
### map options
@click.option(
    "--enzymes",
    type=str,
    default="DpnII,HinfI",
    help='comma-separated list of Hi-C restriction enzymes (dflt: DpnII,HinfI). The default enzymes where chosen in regard of the Arima Hi-C kit.'
)
@click.option(
    "--binnings",
    type=INT_LIST,
    default="10000",
    help='comma separated bin sizes in bp in which each map is generated (dflt: 10000).'
)
@click.option(
    "--format",
    type=click.Choice(["cool", "mcool"]),
    default="cool",
    help='output format of generated maps in each binning: [cool/mcool] (dflt: cool).'
)
### describe options
@click.option(
    "--chroms_list",
    default = "",
    type=str,
    help='comma-separated chromosome list. Uncompatible with the --chroms parameter, which will take priority for providing sequence type information.'
)
@click.option(
    "--mitochondria",
    type=str,
    help="mitochondria (or other organelle) reference. If provided, the organelle will be considered in the signals display. Uncompatible with the --chroms parameter, which will take priority for providing sequence type information."
)
@click.option(
    "--min_chrom_size",
    default = 100000,
    type=int,
    help="Chromosome minimal size for display."
)
@click.option(
    "--chromstart",
    type=str,
    default = "NC_",
    help='Three first character of the chromosome identifier (dftl: "NC_"). All the sequences which id start with chromstart will be considered as chromosomes.'
)
@click.option(
    "--mcool",
    type=str,
    help="path to mcool file of contigs to classify. The program will compute and classify the contact profiles of contigs not referenced in the chromosome info file. Requires --binning."
)
@click.option(
    "--cool",
    type=str,
    default = "",
    help="path to cool file of contigs to classify. The program will compute the descriptive statistics of all contigs not referenced in the chromosome info file if the desired contigs are not provided via the --contigs option."
)
@click.option(
    "--description_binning",
    type=int,
    default="10000",
    help='Binning size for descriptive statistic of HiC generation (dflt: 10000).'
)
@click.option(
    "--formats",
    type=STR_LIST,
    default="pdf",
    help='comma-separated format list for output (image format required).',
)
### processing options
@click.option(
    "-t",
    "--threads",
    type=int,
    default=8,
    help='number of threads to launch for each subtask (dflt: 8)'
)
@click.option(
    "--no_tmp",
    is_flag=True,
    default=False,
    help="if selected, all the temporary files will be discarded. (dflt: False)"
)
@click.option(
    "--sbatch",
    is_flag=True,
    default=False,
    help="if selected, all the bash script will be launched as individual jobs on a SLURM distribution."
)
@click.option(
    "--sbtach_partition",
    default='dedicated',
    type=str,
    help="partition requested for sbatch."
)
@click.option(
    "--sbtach_qos",
    default= 'fast',
    type=str,
    help="quality of service required for sbatch."
)
@click.option(
    "--sbtach_mem",
    default='40G',
    type=str,
    help="minimum amount of real memory requested for sbatch."
)
@click.option(
    "--sbatch_ncpus",
    default=30,
    type=int,
    help="number of cpus required per task fro sbatch."
)
def pipeline(**args):
    kpipeline.pipeline(**args)