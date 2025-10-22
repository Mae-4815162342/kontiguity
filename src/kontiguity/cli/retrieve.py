import click

from .types import *

import kontiguity.retrieve as kretrieve

@click.command("retrieve")
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
    "-i",
    "--index",
    type=str,
    help="path to the reference genome index."
)
@click.option(
    "--min-size",
    type=int,
    default=1000,
    help='minimum size of the kept contigs in bp.'
)
@click.option(
    "--wgs",
    type=PAIR_LIST,
    help="path to the WGS fastq(s). If paired, provide both fastqs comma-separated."
)
@click.option(
    "--table",
    type=str,
    help='path to a csv table providing the data parameters (Mandatory column heads: ["name", "index", "wgs"]).'
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
def retrieve(**args):
    kretrieve.retrieve(**args)