import click

import kontiguity.pipeline as kpipeline

@click.command("pipeline")
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
    type=str,
    help="path to the WGS fastq(s) OR SRA accession. If paired and local, provide both fastqs comma-separated."
)
@click.option(
    "--hic",
    type=str,
    help="path to the Hi-C fastq(s) OR SRA accession. If paired and local, provide both fastqs comma-separated."
)
@click.option(
    "--table",
    type=str,
    help='path to a csv table providing the data parameters (Mandatory column heads: ["name", "ref", "wgs", "hic"]).'
)
@click.option(
    "--start",
    type=str,
    help='starting point of the pipeline. If some processes have been started at this step, the pipeline will pick-up where it left, except in the --no-pickup mode. Requires the arboresence of any previous step to exist. To select in ["load", "retrieve", "map", "classify"] (dflt: load).'
)
@click.option(
    "--stop",
    type=str,
    help='stopping point of the pipeline. To select in ["load", "retrieve", "map", "classify"] (dflt: classify).'
)
@click.option(
    "--no-pickup",
    is_flag=True,
    default=False,
    help="if selected, will restart avery process of the selected starting point in the pipeline."
)
def pipeline(**args):
    click.echo(args)
    kpipeline.pipeline(**args)