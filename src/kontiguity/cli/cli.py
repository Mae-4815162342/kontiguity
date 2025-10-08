import click

from . import (
    load,
    retrieve,
    map,
    classify,
    pipeline
)

@click.group()
@click.pass_context
def cli(ctx):
    pass

cli.add_command(load.load)
cli.add_command(retrieve.retrieve)
cli.add_command(map.map)
cli.add_command(classify.classify)
cli.add_command(pipeline.pipeline)

if __name__ == "__main__":
    cli()