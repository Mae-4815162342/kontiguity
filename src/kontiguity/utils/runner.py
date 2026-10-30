"""Paths, external commands and step bookkeeping shared by every step."""
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import logging
import time


class StepError(RuntimeError):
    pass


# bulky intermediate files (reads, read bins, assemblies, Hi-C alignments, contact maps) can live in a
# separate space (--tmp-dir); registered per output directory, set from the configuration
_TMP = {}


def set_tmp_dir(outdir, tmp_dir):
    """Registers the temporary space of an output directory ('$VAR' and '~' are expanded at run time,
    so '$SCRATCH/kontiguity' works inside SLURM jobs). None: intermediates stay in the output tree."""
    outdir = os.path.abspath(outdir)
    if not tmp_dir:
        _TMP.pop(outdir, None)
        return None
    root = os.path.abspath(os.path.expanduser(os.path.expandvars(str(tmp_dir))))
    if "$" in root:
        raise StepError(f"--tmp-dir {tmp_dir}: environment variable not defined here")
    tag = hashlib.sha1(outdir.encode()).hexdigest()[:6]
    _TMP[outdir] = os.path.join(root, f"kontiguity_{os.path.basename(outdir)}_{tag}")
    return _TMP[outdir]


def tmp_root(outdir):
    return _TMP.get(os.path.abspath(outdir))


class Layout:
    """Output tree of one unit (species x reference).

    With --tmp-dir, reads/, logan/, the retrieve working files, the hicstuff working files and the
    .mcool go to <tmp-dir>/kontiguity_<outdir name>_<hash>/<species>/[<reference>/] instead
    (hic/<unit>.mcool is then a link to the map stored there).

    <outdir>/
      samples.tsv, kontiguity_config.yaml, logs/, slurm/
      <species>/
        reads/wgs/  reads/hic/         downloaded reads (shared by every reference)
        logan/                         Logan contigs, one file per WGS run
        .status/load.done
        <reference>/
          reference/   ref.fa, chromosomes.tsv
          retrieve/    unmapped reads, de novo assemblies, contigs_raw.fa, contigs_raw.tsv
          contigs/     contigs.fa, contigs_info.tsv, genome.fa (= new reference), genome_chromosomes.tsv
          hic/         <unit>.mcool, hicstuff logs
          classification/  contigs_report.tsv, contigs/<id>.json, plots/
          .status/<step>.done
    """

    def __init__(self, outdir, row):
        self.outdir = os.path.abspath(outdir)
        self.row = row
        self.species = row["species"]
        self.unit = row["unit"]
        self.species_dir = os.path.join(self.outdir, row["species"])
        self.unit_dir = os.path.join(self.species_dir, row["ref_slug"])
        # bulky files: in the temporary space when --tmp-dir is given, otherwise in the output tree
        self.tmp = tmp_root(self.outdir)
        self.species_tmp = os.path.join(self.tmp, row["species"]) if self.tmp else self.species_dir
        self.unit_tmp = os.path.join(self.species_tmp, row["ref_slug"]) if self.tmp else self.unit_dir
        self.wgs_dir = os.path.join(self.species_tmp, "reads", "wgs")
        self.hic_dir_reads = os.path.join(self.species_tmp, "reads", "hic")
        self.logan_dir = os.path.join(self.species_tmp, "logan")
        self.ref_dir = os.path.join(self.unit_dir, "reference")
        self.retrieve_dir = os.path.join(self.unit_dir, "retrieve")
        self.contigs_dir = os.path.join(self.unit_dir, "contigs")
        self.hic_dir = os.path.join(self.unit_dir, "hic")
        self.class_dir = os.path.join(self.unit_dir, "classification")
        # read bins, indexes, assemblies; hicstuff alignments and pairs
        self.retrieve_work = os.path.join(self.unit_tmp, "retrieve") if self.tmp else self.retrieve_dir
        self.hic_work = os.path.join(self.unit_tmp, "hic", "hicstuff") if self.tmp else os.path.join(self.hic_dir, "hicstuff")
        self.filter_work = os.path.join(self.unit_tmp, "filter") if self.tmp else os.path.join(self.contigs_dir, "tmp")

    # files -------------------------------------------------------------
    @property
    def ref_fa(self): return os.path.join(self.ref_dir, "ref.fa")
    @property
    def ref_chroms(self): return os.path.join(self.ref_dir, "chromosomes.tsv")
    @property
    def raw_contigs(self): return os.path.join(self.retrieve_dir, "contigs_raw.fa")
    @property
    def raw_info(self): return os.path.join(self.retrieve_dir, "contigs_raw.tsv")
    @property
    def contigs_fa(self): return os.path.join(self.contigs_dir, "contigs.fa")
    @property
    def contigs_info(self): return os.path.join(self.contigs_dir, "contigs_info.tsv")
    @property
    def genome_fa(self): return os.path.join(self.contigs_dir, "genome.fa")
    @property
    def genome_chroms(self): return os.path.join(self.contigs_dir, "genome_chromosomes.tsv")
    @property
    def mcool(self):
        """The contact map as seen by classify/describe (a link to mcool_store when --tmp-dir is used)."""
        return os.path.join(self.hic_dir, f"{self.row['species']}.{self.row['ref_slug']}.mcool")
    @property
    def mcool_store(self):
        """Where 'map' writes the contact map."""
        return os.path.join(os.path.join(self.unit_tmp, "hic") if self.tmp else self.hic_dir,
                            f"{self.row['species']}.{self.row['ref_slug']}.mcool")
    @property
    def report(self): return os.path.join(self.class_dir, "contigs_report.tsv")

    def step_dir(self, step):
        return self.species_dir if step == "load" else self.unit_dir

    def marker(self, step, state="done"):
        return os.path.join(self.step_dir(step), ".status", f"{step}.{state}")

    def is_done(self, step):
        return os.path.isfile(self.marker(step))

    def mark(self, step, state="done", info=None):
        d = os.path.dirname(self.marker(step))
        os.makedirs(d, exist_ok=True)
        for s in ("done", "failed", "running"):
            if s != state and os.path.isfile(self.marker(step, s)):
                os.remove(self.marker(step, s))
        with open(self.marker(step, state), "w") as fh:
            json.dump(dict(time=datetime.datetime.now().isoformat(timespec="seconds"), **(info or {})), fh, indent=1)

    def log(self, step):
        d = os.path.join(self.step_dir(step), "logs")
        os.makedirs(d, exist_ok=True)
        return os.path.join(d, f"{step}.log")


def run(cmd, log_path, logger=None, cwd=None, env=None, check=True):
    """Runs a shell command (bash, pipefail), appending its output to log_path."""
    logger = logger or logging.getLogger("kontiguity")
    os.makedirs(os.path.dirname(log_path) or ".", exist_ok=True)
    logger.info(f"$ {cmd}")
    t0 = time.time()
    with open(log_path, "a") as log:
        log.write(f"\n# {datetime.datetime.now():%Y-%m-%d %H:%M:%S}  $ {cmd}\n")
        log.flush()
        p = subprocess.run(["bash", "-o", "pipefail", "-c", cmd], stdout=log, stderr=subprocess.STDOUT,
                           cwd=cwd, env=env)
    if p.returncode != 0 and check:
        tail = ""
        try:
            with open(log_path) as fh:
                tail = "".join(fh.readlines()[-15:])
        except OSError:
            pass
        raise StepError(f"command failed (exit {p.returncode}) after {time.time() - t0:.0f}s: {cmd}\n"
                        f"--- last lines of {log_path} ---\n{tail}")
    return p.returncode


def require(*tools):
    missing = [t for t in tools if shutil.which(t) is None]
    if missing:
        raise StepError(f"missing executable(s) in PATH: {', '.join(missing)} "
                        f"(install them with the conda environment, see README)")


def gz_writer_cmd(threads=4):
    return f"pigz -p {threads}" if shutil.which("pigz") else "gzip"
