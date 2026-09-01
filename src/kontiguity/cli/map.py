import click

from .types import *

import kontiguity.map as kmap

@click.command("map")
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
    "--hic",
    type=FASTQ,
    help="path to the Hi-C fastq(s). If paired, provide both fastqs comma-separated."
)
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
    "--table",
    type=str,
    help='path to a csv table providing the data parameters (Mandatory column heads: ["name", "index", "hic", "enzymes", "binnings"]).'
)
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
    help="number of cpus required per task from sbatch."
)
def map(**args):
    kmap.map(**args)