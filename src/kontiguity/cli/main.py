"""kontiguity command-line interface.

    kontiguity dtol       build a sample table from the Darwin Tree of Life portal
    kontiguity load       references, Logan contigs, reads           (per species)
    kontiguity retrieve   Logan contigs or read bin -> de novo contigs (per species x reference)
    kontiguity filter     re-alignment, sorting, new reference
    kontiguity map        Hi-C contact map on the new reference (hicstuff + cooler)
    kontiguity classify   contact-profile classification of the contigs
    kontiguity run        all steps at once (local or SLURM)
    kontiguity status     progress of every unit
"""
import copy
import os
import sys

import click

from kontiguity import __version__
from kontiguity.settings import STEPS, load_config, save_config
from kontiguity.samples import resolve_samples, load_samples, select_units
from kontiguity.utils.logging_setup import get_logger


class CSV(click.ParamType):
    name = "a,b,c"

    def __init__(self, cast=str):
        self.cast = cast

    def convert(self, value, param, ctx):
        if isinstance(value, (list, tuple)):
            return list(value)
        try:
            return [self.cast(v) for v in str(value).split(",") if v != ""]
        except ValueError:
            self.fail(f"{value!r} is not a comma-separated list of {self.cast.__name__}", param, ctx)


def _opts(*decorators):
    def deco(f):
        for d in reversed(decorators):
            f = d(f)
        return f
    return deco


# ------------------------------------------------------------------ shared option groups
input_options = _opts(
    click.option("-o", "--outdir", required=True, type=click.Path(file_okay=False), help="Output directory (created if needed)."),
    click.option("--table", type=click.Path(exists=True, dir_okay=False),
                 help="Sample table (CSV/TSV) with Darwin Tree of Life headers: name, ref, wgs, hic "
                      "[+ layout_WGS, protocol_WGS, layout_Hi-C, protocol_Hi-C, chroms, enzymes, wgs_tech]. See README."),
    click.option("-n", "--name", help="Species / experiment name (direct input)."),
    click.option("-r", "--ref", help="Reference(s): assembly accession(s) (GCA_...) or FASTA path(s), comma-separated."),
    click.option("--wgs", help="WGS ('assembly') runs: SRA/ENA run accessions or local FASTQ, comma-separated; "
                               "paired local files joined by ':' (R1.fq.gz:R2.fq.gz)."),
    click.option("--hic", help="Hi-C runs: run accessions or local paired FASTQ 'R1.fq.gz:R2.fq.gz', comma-separated."),
    click.option("--wgs-tech", help="Technology of the --wgs runs (comma-separated, one per run or one for all): "
                                    "hifi | clr | ont | short | auto [auto]."),
    click.option("--hic-protocol", help="Hi-C protocol of the --hic runs (e.g. 'Arima v2', 'Omni-C'), used to choose the enzymes."),
    click.option("--chroms", type=click.Path(exists=True, dir_okay=False),
                 help="Sequence table of a local reference (CSV/TSV: id, sequence_type[, sequence_name]); "
                      "sequence_type: chromosome | organelle | unplaced | ..."),
    click.option("-u", "--unit", "units", multiple=True,
                 help="Restrict to these units ('Species_name/REF'), species or names (repeatable)."),
)
common_options = _opts(
    click.option("-t", "--threads", type=int, help="Threads per job [8]."),
    click.option("--config", "user_config", type=click.Path(exists=True, dir_okay=False),
                 help="YAML file overriding any default parameter (see README, 'Configuration')."),
    click.option("--keep-tmp", is_flag=True, default=None, help="Keep intermediate files (unmapped reads, assemblies, BAM...)."),
    click.option("--tmp-dir", help="Space for the bulky files (reads, Logan contigs, read bins, assemblies, Hi-C "
                                    "alignments, .mcool), e.g. '$SCRATCH/kontiguity'; must be visible from every node."),
    click.option("--no-tmp", is_flag=True, default=None,
                 help="Delete the reads, read bins, Hi-C alignments and the .mcool as soon as no later step needs them."),
    click.option("--force", is_flag=True, help="Redo steps already marked as done."),
    click.option("-v", "--verbose", is_flag=True, help="Also print the log on the terminal."),
)
logan_options = _opts(
    click.option("--no-logan", is_flag=True, default=None, help="Do not use Logan: always download reads and assemble."),
    click.option("--logan-keep-reads", is_flag=True, default=None,
                 help="Download and assemble the reads even when Logan contigs exist (both are used)."),
    click.option("--logan-url", help="Logan contigs URL template [{acc} placeholder]."),
)
load_options = _opts(
    click.option("--download-source", type=click.Choice(["auto", "sra", "ena"]),
                 help="Reads: 'sra' = wget from the NCBI SRA open-data mirror + local fasterq-dump; "
                      "'ena' = ENA fastq.gz; 'auto' = sra then ena [auto]."),
    click.option("--keep-sra", is_flag=True, default=None, help="Keep the downloaded .sra files."),
    click.option("--parallel-downloads", type=int, help="Simultaneous downloads [2]."),
)
retrieve_options = _opts(
    click.option("--aligner", type=click.Choice(["minimap2", "bowtie2"]), help="Read aligner for the read bin [minimap2]."),
    click.option("--min-contig-len", type=int, help="Minimal contig length in bp [1000]."),
    click.option("--megahit-opts", help="Extra MEGAHIT options (short reads), quoted."),
    click.option("--flye-opts", help="Extra Flye options (long reads), quoted."),
)
filter_options = _opts(
    click.option("--contigs", type=click.Path(exists=True, dir_okay=False),
                 help="Use your own contigs (FASTA) as candidates instead of the 'retrieve' output."),
    click.option("--max-ref-cov", type=float, help="Keep contigs with less than this fraction aligned to the reference [0.5]."),
    click.option("--min-identity", type=float, help="Minimal identity of a reference alignment to count [0.90]."),
    click.option("--max-contig-len", type=int, help="Maximal contig length in bp [2000000]."),
    click.option("--no-dedup", is_flag=True, default=None, help="Do not remove redundant contigs."),
    click.option("--minimap2-preset", help="minimap2 preset for contig re-alignment [asm20]."),
)
map_options = _opts(
    click.option("--genome", type=click.Path(exists=True, dir_okay=False),
                 help="Use your own (augmented) reference FASTA instead of the 'filter' output; sequence types from "
                      "--chroms, otherwise sequences >= host-min-len are chromosomes."),
    click.option("--enzymes", help="Hi-C enzymes for hicstuff ('DpnII,HinfI', 'dnase', chunk size...) [auto: from protocol]."),
    click.option("--hic-aligner", type=click.Choice(["bowtie2", "minimap2", "bwa"]), help="hicstuff aligner [bowtie2]."),
    click.option("--hic-mapping", type=click.Choice(["normal", "iterative", "cutsite"]), help="hicstuff mapping mode [cutsite]."),
    click.option("--hic-filter", is_flag=True, default=None, help="hicstuff --filter (removes loops/uncuts; kills intra contacts of small contigs)."),
    click.option("--no-duplicates", is_flag=True, default=None, help="Keep PCR duplicates."),
    click.option("--min-mapq", type=int, help="Minimal mapping quality [30]."),
    click.option("--resolutions", type=CSV(int), help="mcool resolutions, the first one is the base [5000,10000,50000,100000]."),
    click.option("--no-balance", is_flag=True, default=None, help="Do not balance the matrices."),
)
classify_options = _opts(
    click.option("--mcool", type=click.Path(exists=True, dir_okay=False),
                 help="Classify the contigs of an existing .mcool (sequence types from --chroms or by length)."),
    click.option("--model", help="Classifier JSON [builtin]."),
    click.option("--threshold", type=float, help="Probability threshold [0.5]."),
    click.option("--min-depth", type=int, help="Minimal trans Hi-C contacts to evaluate a contig [500]."),
    click.option("--min-effect", type=float, help="Minimal signal amplitude of a strong candidate [0.30]."),
    click.option("--resolution", type=int, help="Resolution used for the profiles [5000, the training resolution]."),
    click.option("--host-min-len", type=int, help="Minimal host chromosome length [100000]."),
    click.option("--recovered-only", is_flag=True, default=None, help="Classify only recovered contigs (not reference scaffolds/organelles)."),
    click.option("--no-contig-files", is_flag=True, default=None, help="Do not write one JSON file per contig."),
    click.option("--plots", type=click.Choice(["none", "candidates", "all"]), help="Contact-profile plots [candidates]."),
)
describe_options = _opts(
    click.option("--known", type=click.Path(exists=True, dir_okay=False),
                 help="FASTA of known sequences (e.g. the 2-micron plasmid) to look for among the recovered contigs."),
    click.option("--describe-contigs", type=click.Choice(["recovered", "candidates", "all"]),
                 help="Contigs drawn by describe [recovered: recovered contigs with Hi-C contacts]."),
    click.option("--describe-max", type=int, help="Maximal number of contigs drawn [50]."),
    click.option("--describe-binnings", type=CSV(int), help="Resolutions of the describe figures [10000]."),
    click.option("--formats", type=CSV(), help="Figure formats [png]."),
    click.option("--no-describe-plots", is_flag=True, default=None, help="Only the tables (recovered contigs, known hits)."),
)
slurm_options = _opts(
    click.option("--slurm/--no-slurm", default=None, help="Submit SLURM jobs instead of running locally (--no-slurm: force local, used inside jobs)."),
    click.option("--slurm-partition"), click.option("--slurm-qos"), click.option("--slurm-account"),
    click.option("--slurm-res", multiple=True,
                 help="Resources of a step: 'step=cpus,mem,time', e.g. 'map=32,96G,72:00:00' (repeatable)."),
    click.option("--slurm-extra", multiple=True, help="Extra #SBATCH line, e.g. '--constraint=avx2' (repeatable)."),
    click.option("--slurm-setup", multiple=True, help="Shell line run first in every job, e.g. 'source activate kontiguity' (repeatable)."),
    click.option("--slurm-wait", is_flag=True, default=None, help="Wait for the submitted jobs to finish."),
    click.option("--dry-run", is_flag=True, help="SLURM: write the job scripts without submitting."),
)

# CLI option -> configuration key(s)
MAP = {
    "threads": [("threads",)], "keep_tmp": [("keep_tmp",)], "tmp_dir": [("tmp_dir",)], "no_tmp": [("no_tmp",)],
    "logan_url": [("logan", "url")],
    "download_source": [("load", "download_source")], "keep_sra": [("load", "keep_sra")],
    "parallel_downloads": [("load", "parallel")],
    "aligner": [("retrieve", "aligner")], "megahit_opts": [("retrieve", "megahit_opts")],
    "flye_opts": [("retrieve", "flye_opts")],
    "min_contig_len": [("retrieve", "min_contig_len"), ("filter", "min_contig_len")],
    "max_ref_cov": [("filter", "max_ref_cov")], "min_identity": [("filter", "min_identity")],
    "max_contig_len": [("filter", "max_contig_len")], "minimap2_preset": [("filter", "minimap2_preset")],
    "enzymes": [("map", "enzymes")], "hic_aligner": [("map", "aligner")], "hic_mapping": [("map", "mapping")],
    "hic_filter": [("map", "filter_events")], "min_mapq": [("map", "min_mapq")], "resolutions": [("map", "resolutions")],
    "model": [("classify", "model")], "threshold": [("classify", "threshold")], "min_depth": [("classify", "min_depth")],
    "min_effect": [("classify", "min_effect")], "resolution": [("classify", "resolution")],
    "host_min_len": [("classify", "host_min_len")], "plots": [("classify", "plots")],
    "known": [("describe", "known")], "describe_contigs": [("describe", "contigs")],
    "describe_max": [("describe", "max_contigs")], "describe_binnings": [("describe", "binnings")],
    "formats": [("describe", "formats")],
    "slurm": [("slurm", "enabled")], "slurm_partition": [("slurm", "partition")], "slurm_qos": [("slurm", "qos")],
    "slurm_account": [("slurm", "account")], "slurm_wait": [("slurm", "wait")],
}
NEGATED = {"no_logan": ("logan", "enabled"), "logan_keep_reads": ("logan", "skip_reads_if_available"),
           "no_dedup": ("filter", "dedup"), "no_duplicates": ("map", "duplicates"), "no_balance": ("map", "balance"),
           "recovered_only": ("classify", "include_reference_scaffolds"), "no_contig_files": ("classify", "contig_files"),
           "no_describe_plots": ("describe", "plots")}


def overrides(kw):
    ov = {}

    def put(path, v):
        d = ov
        for k in path[:-1]:
            d = d.setdefault(k, {})
        d[path[-1]] = v

    for k, v in kw.items():
        if v is None or v == () or v is False and k not in NEGATED:
            continue
        if k in MAP:
            for p in MAP[k]:
                put(p, v)
        elif k in NEGATED and v:
            put(NEGATED[k], False)
    for spec in kw.get("slurm_res") or ():
        step, _, vals = spec.partition("=")
        if step not in STEPS:
            raise click.BadParameter(f"--slurm-res: unknown step '{step}'")
        cpus, mem, time = (vals.split(",") + [None, None, None])[:3]
        put(("slurm", "resources", step), {k: v for k, v in dict(cpus=int(cpus) if cpus else None, mem=mem, time=time).items() if v})
    if kw.get("slurm_extra"):
        put(("slurm", "extra"), list(kw["slurm_extra"]))
    if kw.get("slurm_setup"):
        put(("slurm", "setup"), list(kw["slurm_setup"]))
    return ov


def prepare(kw):
    try:
        return _prepare(kw)
    except (FileNotFoundError, ValueError) as e:
        raise click.ClickException(str(e))


def _prepare(kw):
    samples = resolve_samples(kw["outdir"], kw.get("table"), kw.get("name"), kw.get("ref"), kw.get("wgs"), kw.get("hic"),
                              kw.get("wgs_tech"), kw.get("chroms"), None, kw.get("hic_protocol"), list(kw.get("units") or []))
    cfg = load_config(kw["outdir"], kw.get("user_config"), overrides(kw))
    if cfg.get("no_tmp") and cfg.get("keep_tmp"):
        raise ValueError("--no-tmp and --keep-tmp are incompatible (check kontiguity_config.yaml if you did not "
                         "give both: options are remembered from previous commands)")
    # the execution mode is decided per invocation and never persisted: SLURM jobs read the saved
    # configuration and must run their unit locally (otherwise they would re-submit themselves)
    saved = copy.deepcopy(cfg)
    saved["slurm"]["enabled"] = False
    saved["slurm"]["wait"] = False
    save_config(kw["outdir"], saved)
    if os.environ.get("SLURM_JOB_ID") and kw.get("units") and not kw.get("slurm"):
        cfg["slurm"]["enabled"] = False
    logger = get_logger(kw["outdir"], verbose=kw.get("verbose", False))
    return samples, cfg, logger


def execute(first, last, kw):
    from kontiguity.pipeline import run_pipeline_local, submit
    samples, cfg, logger = prepare(kw)
    logger.info(f"kontiguity {__version__}: steps {first}->{last} on {len(samples)} unit(s): {', '.join(samples.unit)}")
    ui = {k: kw.get(k) for k in ("contigs", "genome", "mcool", "chroms") if kw.get(k)}
    if cfg["slurm"]["enabled"]:
        failed = submit(samples, kw["outdir"], cfg, logger, first, last, kw.get("force", False), kw.get("dry_run", False), ui)
        click.echo(f"jobs written to {kw['outdir']}/slurm/" + (" (dry run, nothing submitted)" if kw.get("dry_run") else ""))
    else:
        failed = run_pipeline_local(samples, kw["outdir"], cfg, logger, first, last, kw.get("force", False), ui)
    if failed:
        click.echo(f"FAILED: {', '.join(failed)} (see {kw['outdir']}/logs/kontiguity.log and <unit>/.status/)", err=True)
        sys.exit(1)


@click.group(context_settings=dict(help_option_names=["-h", "--help"], max_content_width=110))
@click.version_option(__version__, prog_name="kontiguity")
def cli():
    """Recovery of extrachromosomal eukaryotic contigs and Hi-C contact classification (plasmid detection)."""


def _step_command(step, *groups, doc=""):
    @cli.command(step, help=doc)
    @input_options
    @common_options
    @slurm_options
    @_opts(*groups)
    def cmd(**kw):
        execute(step, step, kw)
    return cmd


_step_command("load", logan_options, load_options,
              doc="Download references (ENA), Logan contigs and reads (wget + local fasterq-dump). One job per species.")
_step_command("retrieve", logan_options, retrieve_options,
              doc="Candidate contigs: Logan contigs when available, otherwise de novo assembly of the reads that do not "
                  "align on the reference (MEGAHIT for short reads, Flye --meta for long reads).")
_step_command("filter", filter_options,
              doc="Re-align the contigs on the reference, keep those that do not match it, remove redundancy and build "
                  "the new reference (reference + kept contigs).")
_step_command("map", map_options,
              doc="Map the Hi-C reads on the new reference with hicstuff and build a multi-resolution .mcool with cooler.")
_step_command("classify", classify_options,
              doc="Classify the contact profiles of the contigs (plasmid-like 'hitchhiking' signal) and write the per-contig report.")
_step_command("describe", describe_options,
              doc="Describe the recovered contigs: table of all recovered contigs, search of known sequences (--known) and "
                  "descriptive Hi-C figures (contact maps, signals, mini matrices).")


@cli.command("run")
@input_options
@common_options
@slurm_options
@_opts(logan_options, load_options, retrieve_options, filter_options, map_options, classify_options, describe_options)
@click.option("--from", "first", type=click.Choice(STEPS), default="load", show_default=True, help="First step.")
@click.option("--to", "last", type=click.Choice(STEPS), default="describe", show_default=True, help="Last step.")
def run_cmd(first, last, **kw):
    """Run the whole pipeline (or the --from/--to range): load -> retrieve -> filter -> map -> classify -> describe.
    With --slurm, every step of every unit is submitted at once, chained by job dependencies."""
    if STEPS.index(first) > STEPS.index(last):
        raise click.BadParameter("--from is after --to")
    execute(first, last, kw)


@cli.command("dtol")
@click.option("-o", "--output", required=True, help="Output table (.tsv or .csv).")
@click.option("--species", multiple=True, help="Keep species matching this regular expression (repeatable).")
@click.option("--ref-policy", type=click.Choice(["all", "latest", "annotated"]), default="all", show_default=True,
              help="Which assemblies of a species become references.")
@click.option("--primary-only", is_flag=True, help="Drop alternate-haplotype assemblies.")
@click.option("--wgs-tech", type=CSV(), help="Keep only WGS runs of these technologies (hifi,clr,ont,short).")
@click.option("--max-records", type=int, help="Only fetch the first N portal records (tests).")
@click.option("--from-json", type=click.Path(exists=True), help="Parse a saved portal dump instead of querying the API.")
@click.option("--save-json", help="Save the raw portal records (to re-use with --from-json).")
@click.option("--url", default=None, help="Portal API URL.")
@click.option("-t", "--threads", type=int, default=8, show_default=True)
@click.option("-v", "--verbose", is_flag=True)
def dtol_cmd(output, species, ref_policy, primary_only, wgs_tech, max_records, from_json, save_json, url, threads, verbose):
    """Build a sample table from the Darwin Tree of Life portal: species with an assembly, Hi-C and WGS runs."""
    import os
    from kontiguity.dtol import dtol_table
    from kontiguity.settings import DEFAULTS
    logger = get_logger(os.path.dirname(os.path.abspath(output)), verbose=verbose)
    df = dtol_table(output, url or DEFAULTS["urls"]["dtol"], ref_policy, primary_only, wgs_tech, list(species), max_records,
                    threads, save_json, from_json, logger)
    click.echo(f"{df.name.nunique()} species, {len(df)} (species, assembly) rows -> {output}")


@cli.command("status")
@click.option("-o", "--outdir", required=True, type=click.Path(exists=True, file_okay=False))
@click.option("-u", "--unit", "units", multiple=True)
def status_cmd(outdir, units):
    """Progress of every unit (done / running / failed / pending) and number of candidates."""
    import pandas as pd
    from kontiguity.pipeline import status_table
    try:
        df = status_table(select_units(load_samples(outdir), list(units)), outdir)
    except (FileNotFoundError, ValueError) as e:
        raise click.ClickException(str(e))
    with pd.option_context("display.width", 200, "display.max_rows", 500):
        click.echo(df.to_string(index=False))




if __name__ == "__main__":
    cli()
