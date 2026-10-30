"""Default parameters of every kontiguity step.

Precedence (lowest -> highest):
    DEFAULTS  <  <outdir>/kontiguity_config.yaml (saved by previous commands)
              <  --config user.yaml
              <  command-line options

The effective configuration is written back to <outdir>/kontiguity_config.yaml so
that every later command (and every SLURM job) uses exactly the same parameters.
"""
import copy
import os

import yaml

STEPS = ["load", "retrieve", "filter", "map", "classify", "describe"]

DEFAULTS = {
    "threads": 8,
    "keep_tmp": False,
    # bulky intermediates (reads, Logan contigs, read bins, assemblies, Hi-C alignments, .mcool) are
    # written to <tmp_dir>/kontiguity_<outdir name>_<hash>/ instead of the output tree; '$VAR' allowed
    "tmp_dir": None,
    # delete them as soon as no later step needs them (reads, read bins, alignments; .mcool after describe)
    "no_tmp": False,
    # ---------------------------------------------------------------- load
    "load": {
        # where reads are downloaded from:
        #   auto = SRA open-data mirror (wget) + local fasterq-dump, falling back to ENA fastq
        #   sra  = SRA open-data mirror only;  ena = ENA fastq.gz only
        "download_source": "auto",
        "keep_sra": False,
        "compress_fastq": True,
        "wget_tries": 5,
        "wget_timeout": 120,
        # ENA sequence types kept in the reference (others are dropped)
        "sequence_types": ["chromosome", "organelle", "unlocalized", "unplaced", "scaffold", "contig"],
    },
    # ---------------------------------------------------------------- Logan
    "logan": {
        "enabled": True,
        # if Logan has contigs for a WGS run, its reads are not downloaded at all
        "skip_reads_if_available": True,
        "url": "https://s3.amazonaws.com/logan-pub/c/{acc}/{acc}.contigs.fa.zst",
    },
    # ---------------------------------------------------------------- retrieve
    "retrieve": {
        "aligner": "minimap2",            # minimap2 | bowtie2 (bowtie2: short reads only)
        "assembler_short": "megahit",     # megahit
        "assembler_long": "flye",         # flye (metagenome mode, reports circular contigs)
        "megahit_opts": "",
        "flye_opts": "--min-overlap 1000",   # Flye picks a min overlap close to the read N90 otherwise: too strict for read bins
        "min_contig_len": 1000,
        "long_read_max_aligned": 0.5,     # long reads with < 50% aligned also go to the read bin
        "read_circles": True,             # HiFi: circles detected as tandem copies inside reads
        "read_circle_min_len": 1000,
    },
    # ---------------------------------------------------------------- filter
    "filter": {
        "min_contig_len": 1000,
        "max_contig_len": 2_000_000,
        # a contig is kept as 'new' if less than max_ref_cov of its length aligns to the
        # reference with at least min_identity
        "max_ref_cov": 0.5,
        "min_identity": 0.90,
        "minimap2_preset": "asm20",
        "dedup": True,
        "dedup_cov": 0.95,
        "dedup_identity": 0.95,
        "circularity_min_overlap": 20,
        "circularity_max_overlap": 300,
        "circularity_max_mismatch": 0.02,
    },
    # ---------------------------------------------------------------- map
    "map": {
        "enzymes": "auto",               # auto = from the Hi-C protocol (Arima -> DpnII,HinfI)
        "aligner": "bowtie2",            # hicstuff aligner: bowtie2 | minimap2 | bwa
        "mapping": "cutsite",            # hicstuff mapping: normal | iterative | cutsite
        "filter_events": False,          # hicstuff --filter (removes intra-contig contacts of small contigs!)
        "duplicates": True,              # hicstuff --duplicates
        "min_mapq": 30,
        "resolutions": [5000, 10000, 50000, 100000],
        "balance": True,
    },
    # ---------------------------------------------------------------- classify
    "classify": {
        "resolution": 5000,
        "model": "builtin",
        "host_min_len": 100_000,
        "host_types": ["chromosome"],
        "min_depth": 500,                # trans contacts below which a profile is not evaluable
        "threshold": 0.5,
        "min_effect": 0.30,              # minimal signal amplitude (weakest known positive: 0.35)
        "min_intra": 0.005,
        "max_intra": 0.5,
        "community_small": 300_000,      # communities smaller than this are 'molecules'
        "community_frac": 0.2,           # >= this share of contacts with other contigs of a big community -> genome piece
        "cis_frac": 0.10,                # >= this share of contacts in one 25-kb host window -> cis-anchored fragment
        "n_null": 8,
        "contig_files": True,            # one JSON file per contig
        "plots": "candidates",           # none | candidates | all
        "include_reference_scaffolds": True,
        "seed": 0,
    },
    # ---------------------------------------------------------------- describe
    "describe": {
        "known": None,                   # FASTA of known sequences (e.g. 2-micron) to look for among the contigs
        "plots": True,                   # descriptive figures (scripts/describe/describe.py)
        "contigs": "recovered",          # recovered (with Hi-C contacts) | candidates | all
        "max_contigs": 50,               # figures are drawn for at most this many contigs
        "binnings": [10000],
        "formats": ["png"],
        "mini": True,                    # mini contact matrices
    },
    # ---------------------------------------------------------------- SLURM
    "slurm": {
        "enabled": False,
        "partition": None,
        "qos": None,
        "account": None,
        "extra": [],                      # extra '#SBATCH ...' lines, e.g. ["--constraint=avx2"]
        "setup": [],                      # shell lines run before the command (module load / conda activate)
        "wait": False,
        "poll_seconds": 60,
        "resources": {                    # cpus, mem, time per step
            "load":     {"cpus": 4,  "mem": "8G",  "time": "24:00:00"},
            "retrieve": {"cpus": 16, "mem": "64G", "time": "48:00:00"},
            "filter":   {"cpus": 8,  "mem": "32G", "time": "12:00:00"},
            "map":      {"cpus": 16, "mem": "48G", "time": "48:00:00"},
            "classify": {"cpus": 4,  "mem": "16G", "time": "06:00:00"},
            "describe": {"cpus": 2,  "mem": "16G", "time": "06:00:00"},
        },
    },
    # ---------------------------------------------------------------- remote services
    "urls": {
        "sra": "https://sra-pub-run-odp.s3.amazonaws.com/sra/{acc}/{acc}",
        "ena_filereport": ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}"
                           "&result=read_run&fields=run_accession,fastq_ftp,library_layout,instrument_platform&format=tsv"),
        "ena_fasta": "https://www.ebi.ac.uk/ena/browser/api/fasta/{acc}?download=true&gzip=true",
        "dtol": "https://portal.darwintreeoflife.org/api/data_portal",
    },
}

CONFIG_NAME = "kontiguity_config.yaml"


def deep_update(base, new):
    for k, v in (new or {}).items():
        if isinstance(v, dict) and isinstance(base.get(k), dict):
            deep_update(base[k], v)
        elif v is not None:
            base[k] = v
    return base


def load_config(outdir, user_config=None, overrides=None):
    cfg = copy.deepcopy(DEFAULTS)
    saved = os.path.join(outdir, CONFIG_NAME)
    if os.path.isfile(saved):
        with open(saved) as fh:
            deep_update(cfg, yaml.safe_load(fh) or {})
    if user_config:
        with open(user_config) as fh:
            deep_update(cfg, yaml.safe_load(fh) or {})
    deep_update(cfg, overrides or {})
    return cfg


def save_config(outdir, cfg):
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, CONFIG_NAME), "w") as fh:
        yaml.safe_dump(cfg, fh, sort_keys=False)
