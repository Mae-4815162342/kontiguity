import click
from .types import *
import kontiguity.load as kload

@click.command("load")
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
    "-r",
    "--ref",
    type=str,
    help="path to the reference genome fasta OR the GCA reference which will automatically be loaded from ENA database."
)
@click.option(
    "--chroms",
    type=str,
    help='path to a chromosome information file detailing the type of each sequence present in the reference (Mandatory column heads: ["id", "sequence_type", "sequence_name"]). "sequence_type" must be in the ENA database format : ["chromosome", "organelle", ...]. Required only for a local fasta, GCA referenced genomes will have the chromosome.tsv generated.'
)
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
    "--table",
    type=str,
    help='path to a csv table providing the data parameters (Mandatory column heads: ["name", "ref", "wgs", "hic"]). See test_data/load/test_dataset.csv for format example.'
)
@click.option(
    "--dtol",
    is_flag=True,
    default=False,
    help="if selected, a data table will be created and loaded from the Darwin Tree of Life project database."
)
@click.option(
    "-t",
    "--threads",
    type=int,
    default=8,
    help='number of threads to launch for each subtask (dflt: 8)'
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
def load(**args):
    kload.load(**args)