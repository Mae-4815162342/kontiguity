"""Step 'load' (one job per species): references, Logan contigs and reads.

* reference: local FASTA (+ optional chromosome table) or an assembly accession (GCA_...) fetched
  from ENA; sequences are renamed to their accession and typed (chromosome, organelle, ...).
* Logan: for each WGS run with an SRA accession, the Logan contigs are fetched when they exist.
  The reads of such runs are then NOT downloaded (logan.skip_reads_if_available).
* reads: wget of the SRA object from the NCBI open-data mirror, then local fasterq-dump
  (prefetch is not used: it does not work on most SLURM compute nodes); ENA fastq.gz as fallback.
"""
import glob
import os
import shutil
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import requests

from kontiguity import logan
from kontiguity.fasta import read_fasta, write_record, parse_ena_header, _open
from kontiguity.samples import reads_of, is_sra, ASSEMBLY_RE
from kontiguity.utils.runner import Layout, run, require, StepError, gz_writer_cmd

MANIFEST = "manifest.tsv"


# ------------------------------------------------------------------ references
def load_reference(lay, cfg, logger):
    row, log = lay.row, lay.log("load")
    os.makedirs(lay.ref_dir, exist_ok=True)
    if os.path.isfile(lay.ref_fa) and os.path.isfile(lay.ref_chroms):
        logger.info(f"[load:{lay.unit}] reference already present")
        return
    ref = row["ref"]
    if not ref:
        raise StepError(f"{lay.unit}: no reference given")
    if os.path.exists(ref):
        src = ref
    else:  # any other value is treated as an ENA/INSDC assembly accession (GCA_..., GCF_...)
        if not ASSEMBLY_RE.match(ref):
            logger.warning(f"[load:{lay.unit}] '{ref}' is not a file and does not look like GCA_/GCF_; trying ENA anyway")
        require("wget")
        src = os.path.join(lay.ref_dir, "download.fa.gz")
        url = cfg["urls"]["ena_fasta"].format(acc=ref)
        run(f"wget -q --tries={cfg['load']['wget_tries']} --timeout={cfg['load']['wget_timeout']} -O {src}.part '{url}' "
            f"&& mv {src}.part {src}", log, logger)

    user_types = read_chroms_table(row.get("chroms", "")) if row.get("chroms") else None
    keep_types = set(cfg["load"]["sequence_types"])
    rows = []
    with open(lay.ref_fa + ".part", "w") as out:
        for header, seq in read_fasta(src):
            cid, stype, sname, desc = parse_ena_header(header)
            if user_types is not None:
                info = user_types.get(cid) or user_types.get(header.split()[0])
                stype, sname = (info if info else ("unknown", cid))
            elif stype == "unknown" and len(seq) >= cfg["classify"]["host_min_len"]:
                stype = "chromosome"   # untyped local FASTA: long sequences are considered chromosomes
            if stype not in keep_types and stype != "unknown":
                continue
            write_record(out, cid, seq.upper())
            rows.append(dict(id=cid, sequence_type=stype, sequence_name=sname, length=len(seq), origin="reference",
                             description=desc))
    os.replace(lay.ref_fa + ".part", lay.ref_fa)
    pd.DataFrame(rows).to_csv(lay.ref_chroms, sep="\t", index=False)
    if src.endswith("download.fa.gz"):
        os.remove(src)
    n_chr = sum(r["sequence_type"] == "chromosome" for r in rows)
    logger.info(f"[load:{lay.unit}] reference: {len(rows)} sequences ({n_chr} chromosomes)")


def read_chroms_table(path):
    """User chromosome table (csv/tsv) with columns id, sequence_type[, sequence_name]."""
    sep = "\t" if path.endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep, dtype=str)
    if not {"id", "sequence_type"} <= set(df.columns):
        raise StepError(f"{path}: chromosome table needs columns id, sequence_type[, sequence_name]")
    return {r.id: (r.sequence_type, getattr(r, "sequence_name", r.id)) for r in df.itertuples()}


# ------------------------------------------------------------------ reads
def _final_names(outdir, acc):
    return dict(pair=[os.path.join(outdir, f"{acc}_1.fastq.gz"), os.path.join(outdir, f"{acc}_2.fastq.gz")],
                single=[os.path.join(outdir, f"{acc}.fastq.gz")])


def existing_reads(outdir, acc):
    n = _final_names(outdir, acc)
    if all(os.path.isfile(f) for f in n["pair"]):
        return n["pair"]
    if os.path.isfile(n["single"][0]):
        return n["single"]
    return None


def download_sra(acc, outdir, cfg, threads, log, logger):
    """wget the SRA object from the open-data mirror, then fasterq-dump it locally."""
    require("wget", "fasterq-dump")
    work = os.path.join(outdir, f"{acc}.sra_tmp")
    os.makedirs(work, exist_ok=True)
    sra = os.path.join(work, acc)   # fasterq-dump names the output after the file name
    url = cfg["urls"]["sra"].format(acc=acc)
    rc = run(f"wget -q -c --tries={cfg['load']['wget_tries']} --timeout={cfg['load']['wget_timeout']} -O {sra} '{url}'",
             log, logger, check=False)
    if rc != 0 or not os.path.isfile(sra) or os.path.getsize(sra) == 0:
        shutil.rmtree(work, ignore_errors=True)
        return None
    run(f"cd {work} && fasterq-dump --split-3 --skip-technical -e {threads} --temp {work} -O {work} ./{acc}", log, logger)
    fq = sorted(glob.glob(os.path.join(work, "*.fastq")))
    pairs = [f for f in fq if f.endswith(("_1.fastq", "_2.fastq"))]
    gz = gz_writer_cmd(threads)
    names = _final_names(outdir, acc)
    if len(pairs) == 2:
        for src, dst in zip(sorted(pairs), names["pair"]):
            run(f"{gz} -c {src} > {dst}.part && mv {dst}.part {dst}", log, logger)
        files = names["pair"]
    elif fq:
        run(f"cat {' '.join(fq)} | {gz} -c > {names['single'][0]}.part && mv {names['single'][0]}.part {names['single'][0]}",
            log, logger)
        files = names["single"]
    else:
        raise StepError(f"fasterq-dump produced no fastq for {acc}")
    if cfg["load"]["keep_sra"]:
        shutil.move(sra, os.path.join(outdir, f"{acc}.sra"))
    shutil.rmtree(work, ignore_errors=True)
    return files


def download_ena(acc, outdir, cfg, log, logger):
    url = cfg["urls"]["ena_filereport"].format(acc=acc)
    try:
        txt = requests.get(url, timeout=cfg["load"]["wget_timeout"]).text
    except requests.RequestException as e:
        logger.warning(f"ENA filereport failed for {acc}: {e}")
        return None
    lines = [l.split("\t") for l in txt.strip().splitlines()]
    if len(lines) < 2:
        return None
    head = lines[0]
    col = head.index("fastq_ftp") if "fastq_ftp" in head else 1
    urls = [u for u in lines[1][col].split(";") if u]
    if not urls:
        return None
    files = []
    for u in urls:
        dst = os.path.join(outdir, os.path.basename(u))
        u = u if u.startswith(("http://", "https://", "ftp://")) else "https://" + u
        run(f"wget -q -c --tries={cfg['load']['wget_tries']} --timeout={cfg['load']['wget_timeout']} -O {dst}.part '{u}' "
            f"&& mv {dst}.part {dst}", log, logger)
        files.append(dst)
    # ENA names: ACC_1.fastq.gz, ACC_2.fastq.gz (+ optional ACC.fastq.gz singletons)
    pair = [f for f in files if f.endswith(("_1.fastq.gz", "_2.fastq.gz"))]
    return sorted(pair) if len(pair) == 2 else files[:1]


def fetch_reads(acc, outdir, cfg, threads, log, logger):
    os.makedirs(outdir, exist_ok=True)
    have = existing_reads(outdir, acc)
    if have:
        return have, "cached"
    source = cfg["load"]["download_source"]
    if source in ("auto", "sra"):
        files = download_sra(acc, outdir, cfg, threads, log, logger)
        if files:
            return files, "sra"
        logger.warning(f"[load] {acc}: not found on the SRA open-data mirror" + (", trying ENA" if source == "auto" else ""))
    if source in ("auto", "ena"):
        files = download_ena(acc, outdir, cfg, log, logger)
        if files:
            return files, "ena"
    raise StepError(f"could not download reads of {acc} (source={source})")


# ------------------------------------------------------------------ species step
def load_species(lays, cfg, logger, need=("ref", "wgs", "hic")):
    """lays: list of Layout of the units of ONE species.  `need` restricts what is prepared
    (a step run on its own only prepares its own inputs, e.g. 'map' needs the Hi-C reads only)."""
    lay0 = lays[0]
    log = lay0.log("load")
    threads = cfg["threads"]
    if "ref" in need:
        for lay in lays:
            load_reference(lay, cfg, logger)

    # union of the runs of the species
    wgs, hic = {}, {}
    for lay in lays:
        for e in reads_of(lay.row, "wgs"):
            wgs.setdefault(e["id"], e)
        for e in reads_of(lay.row, "hic"):
            hic.setdefault(e["id"], e)

    if "wgs" not in need:
        wgs = {}
    if "hic" not in need:
        hic = {}
    path = os.path.join(lay0.species_dir, "reads", MANIFEST)
    previous = {}
    if os.path.isfile(path):   # keep what a previous (partial) load already prepared
        for m in pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False).to_dict("records"):
            # downloads deleted since (--no-tmp) are prepared again
            if m["source"] != "pending" and not _entry_missing(m):
                previous[(m["kind"], m["run"])] = m
    manifest = []
    # Logan first: runs it covers need no read download
    lg = cfg["logan"]
    for rid, e in wgs.items():
        if ("wgs", rid) in previous:
            manifest.append(previous.pop(("wgs", rid)))
            continue
        if e["files"]:
            manifest.append(dict(kind="wgs", run=rid, files=":".join(e["files"]), layout=e["layout"], tech=e["tech"],
                                 source="local", logan=""))
            continue
        logan_fa, status = None, "disabled"
        if lg["enabled"] and is_sra(rid):
            ok = logan.available(rid, lg["url"])
            status = {True: "available", False: "absent", None: "unknown"}[ok]
            if ok is not False:
                logan_fa = logan.download(rid, lg["url"], lay0.logan_dir, log, logger,
                                          cfg["load"]["wget_tries"], cfg["load"]["wget_timeout"])
                status = "downloaded" if logan_fa else "absent"
            logger.info(f"[load:{lay0.species}] Logan {rid}: {status}")
        need_reads = not (logan_fa and lg["skip_reads_if_available"])
        manifest.append(dict(kind="wgs", run=rid, files="", layout=e["layout"], tech=e["tech"],
                             source="pending" if need_reads else "logan_only", logan=logan_fa or status))

    jobs = [m for m in manifest if m["source"] == "pending"]
    for rid, e in hic.items():
        if ("hic", rid) in previous:
            manifest.append(previous.pop(("hic", rid)))
            continue
        m = dict(kind="hic", run=rid, files=":".join(e["files"]), layout=e["layout"] or "PAIRED", tech="hic",
                 source="local" if e["files"] else "pending", logan="")
        manifest.append(m)
        if not e["files"]:
            jobs.append(m)

    def _dl(m):
        outdir = lay0.wgs_dir if m["kind"] == "wgs" else lay0.hic_dir_reads
        files, src = fetch_reads(m["run"], outdir, cfg, max(1, threads // max(1, cfg["load"].get("parallel", 2))), log, logger)
        m["files"], m["source"] = ":".join(files), src
        m["layout"] = "PAIRED" if len(files) == 2 else "SINGLE"
        logger.info(f"[load:{lay0.species}] {m['kind']} {m['run']}: {len(files)} file(s) from {src}")

    if jobs:
        with ThreadPoolExecutor(max_workers=max(1, cfg["load"].get("parallel", 2))) as ex:
            list(ex.map(_dl, jobs))
    manifest += list(previous.values())
    os.makedirs(os.path.dirname(path), exist_ok=True)
    pd.DataFrame(manifest, columns=["kind", "run", "files", "layout", "tech", "source", "logan"]).to_csv(path, sep="\t", index=False)
    return path


def _entry_missing(m):
    """True if a manifest entry points to downloaded files (reads or Logan contigs) that no longer exist."""
    if m["source"] not in ("local", "logan_only") and any(not os.path.isfile(f) for f in m["files"].split(":") if f):
        return True
    lg = m.get("logan") or ""
    return lg.startswith("/") and not os.path.isfile(lg)


def missing_inputs(lay, kinds):
    """Kinds of downloaded inputs ('wgs', 'hic') of a species that are absent (never loaded or deleted)."""
    path = os.path.join(lay.species_dir, "reads", MANIFEST)
    if not os.path.isfile(path):
        return set(kinds)
    man = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    return {k for k in kinds if k in ("wgs", "hic") and
            any(m["source"] == "pending" or _entry_missing(m) for m in man[man.kind == k].to_dict("records"))}


def read_manifest(lay):
    path = os.path.join(lay.species_dir, "reads", MANIFEST)
    if not os.path.isfile(path):
        raise StepError(f"{path} missing: run 'kontiguity load' first")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
