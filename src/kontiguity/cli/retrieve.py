import click

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
    type=str,
    help="path to the WGS fastq(s). If paired, provide both fastqs comma-separated."
)
@click.option(
    "--table",
    type=str,
    help='path to a csv table providing the data parameters (Mandatory column heads: ["name", "index", "wgs"]).'
)
def retrieve(**args):
    click.echo(args)
    kretrieve.retrieve(**args)