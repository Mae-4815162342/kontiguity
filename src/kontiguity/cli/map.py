import click

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
    type=str,
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
    type=str,
    default="10000",
    help='comma separated bin sizes in bp in which each map is generated (dflt: 10000).'
)
@click.option(
    "--zoomify",
    is_flag=True,
    default=False,
    help="if provided will produce a mcool file instead of separated cools. In such case, the smalest binning in the binnings list must be a common divider of the other values (e.g. 1000,2000,5000)."
)
@click.option(
    "--table",
    type=str,
    help='path to a csv table providing the data parameters (Mandatory column heads: ["name", "index", "hic", "enzymes", "binnings"]).'
)
def map(**args):
    click.echo(args)
    kmap.map(**args)