import click
from .types import *
import kontiguity.describe as kdescribe

@click.command("describe")
@click.option(
    '-n',
    '--name', 
    type=str,
    help="name of the experiment (recommanded: species name. info: spaces are not allowed and will be replaced by _.)"
)
@click.option(
    '-o',
    '--outpath',
    type=str,
    help="output folder path, created if non-existent"
)
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
    "--chroms",
    type=str,
    default="",
    help='path to a chromosome information file detailing the type of each sequence present in the reference (Mandatory column heads: ["id", "sequence_type", "sequence_name"]). "sequence_type" must be in the ENA database format : ["chromosome", "organelle", ...]. Required only for a local fasta, GCA referenced genomes will have the chromosome.tsv generated.'
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
    "--contigs",
    default = "",
    type=str,
    help="comma-separated list of the contigs to describe. If not provided, will evaluate all the contigs that are not identified as chromosomes with a total trans-coverage of at least 1."
)
@click.option(
    "-i",
    "--index",
    type=str,
    help="path to the reference genome index."
)
@click.option(
    "--mcool",
    type=str,
    help="path to mcool file of contigs to classify. The program will compute and classify the contact profiles of contigs not referenced in the chromosome info file. Requires --binning."
)
@click.option(
    "--cool",
    type=str,
    help="path to cool file of contigs to classify. The program will compute the descriptive statistics of all contigs not referenced in the chromosome info file if the desired contigs are not provided via the --contigs option."
)
@click.option(
    "--binnings",
    type=INT_LIST,
    default="10000",
    help='comma separated bin sizes in bp in which each map is generated (dflt: 10000).'
)
@click.option(
    "--formats",
    type=STR_LIST,
    default="pdf",
    help='comma-separated format list for output (image format required).',
)
@click.option(
    "--table",
    type=str,
    help='path to a csv table providing the data parameters (Mandatory column heads: ["name", "ref", "wgs", "hic"]).'
)
@click.option(
    "--mini_only",
    is_flag=True,
    default=False,
    help="if selected, the only output of describe is mini-matrices. For each cool in the dataset, mini-matrices, unbalanced and normalized, will be produce, containing all the chromosomes and contigs with a minimal trans contact signal of 0."
)
@click.option(
    "--no_mini",
    is_flag=True,
    default=False,
    help="if selected, mini matrices single figure are not outputed. Mini matrices are still present in the summary figures."
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
    "-v",
    "--verbose",
    is_flag=True,
    default=False,
    help="if selected, also prints log messages to the terminal in addition to the run's log file (which is always written). Off by default to avoid crowding the terminal on large datasets."
)
@click.option(
    "--sbatch_ncpus",
    default=30,
    type=int,
    help="number of cpus required per task fro sbatch."
)
def describe(**args):
    kdescribe.describe(**args)