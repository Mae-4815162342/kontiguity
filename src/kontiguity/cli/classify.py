import click

import kontiguity.classify as kclassify

@click.command("classify")
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
    "--chroms",
    type=str,
    help='path to a chromosome information file detailing the type of each sequence present in the reference (Mandatory column heads: ["id", "sequence_type", "sequence_name"]). "sequence_type" must be in the ENA database format : ["chromosome", "organelle", ...]. Required only for a local fasta, GCA referenced genomes will have the chromosome.tsv generated.'
)
@click.option(
    "--mcool",
    type=str,
    help="path to mcool file of contigs to classify. The program will compute and classify the contact profiles of contigs not referenced in the chromosome info file. Requires --binning."
)
@click.option(
    "--binnings",
    type=str,
    default="10000",
    help='comma separated bin sizes in bp in which each map is generated (dflt: 10000).'
)
@click.option(
    "--cool",
    type=str,
    help="ppath to cool file of contigs to classify. The program will compute and classify the contact profiles of contigs not referenced in the chromosome info file."
)
@click.option(
    "--model",
    type=str,
    default="models/model.hdf5",
    help="path to a classifier in hdf5 format (dflt: provided model)."
)
@click.option(
    "--param_file",
    type=str,
    default="models/model.json",
    help="path to a json file containing the data preprodessing and model parameters (dflt: model params file)"
)
@click.option(
    "--table",
    type=str,
    help='path to a csv table providing the data parameters (Mandatory column heads: ["name", "ref", "wgs", "hic"]).'
)
def classify(**args):
    click.echo(args)
    kclassify.classify(**args)