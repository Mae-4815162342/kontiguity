"""Sample tables.

Every command accepts its input in one of three ways:
  1. direct identifiers:  --name, --ref, --wgs, --hic  (accessions or local files)
  2. a table (--table) using the Darwin Tree of Life headers
     (name, ref, wgs, hic, layout_WGS, protocol_WGS, layout_Hi-C, protocol_Hi-C, ...)
  3. nothing: the normalized table <outdir>/samples.tsv written by a previous command is reused.

Internally the analysis unit is (species, reference assembly).  Reads belong to the species
and are shared by all its references; Logan contigs are fetched per WGS run.
"""
import os
import re

import pandas as pd

SRA_RE = re.compile(r"^[SED]RR\d+$")
ASSEMBLY_RE = re.compile(r"^GC[AF]_\d+(\.\d+)?$")
SAMPLES_NAME = "samples.tsv"

SAMPLE_COLUMNS = ["unit", "name", "species", "ref", "ref_slug", "chroms",
                  "wgs", "wgs_files", "wgs_layout", "wgs_tech", "wgs_protocol",
                  "hic", "hic_files", "hic_layout", "hic_protocol", "enzymes"]


def slug(text):
    s = re.sub(r"[^A-Za-z0-9.]+", "_", str(text).strip()).strip("_")
    return s or "unnamed"


def is_sra(x):
    return bool(SRA_RE.match(str(x).strip()))


def ref_slug(ref):
    ref = str(ref).strip()
    if os.path.exists(ref):
        b = os.path.basename(ref)
        for ext in (".gz", ".fasta", ".fna", ".fa", ".fas"):
            if b.endswith(ext):
                b = b[: -len(ext)]
        return slug(b)
    return slug(ref)


# ------------------------------------------------------------------ inference helpers
def infer_tech(protocol="", platform="", layout=""):
    """Sequencing technology of a WGS run: hifi | clr | ont | short | auto."""
    p = f"{protocol} {platform}".lower()
    if any(k in p for k in ("hifi", "ccs", "revio", "sequel ii")):
        return "hifi"
    if any(k in p for k in ("nanopore", "oxford", " ont", "minion", "promethion", "gridion")) or p.startswith("ont"):
        return "ont"
    if "clr" in p or ("pacbio" in p and "hifi" not in p):
        return "clr"
    if any(k in p for k in ("illumina", "chromium", "10x", "truseq", "nextera", "novaseq", "hiseq", "nextseq", "bgiseq", "dnbseq", "pcr")):
        return "short"
    if str(layout).upper() == "PAIRED":
        return "short"
    return "auto"


def infer_enzymes(protocol=""):
    """hicstuff --enzyme value from a Hi-C library protocol string."""
    p = str(protocol).lower()
    if "arima" in p:
        return "DpnII,HinfI"
    if any(k in p for k in ("omni-c", "omnic", "omni c", "dovetail", "dnase")):
        return "dnase"
    if "micro-c" in p or "mnase" in p:
        return "mnase"
    if any(k in p for k in ("phase", "proximo", "sau3ai", "dpnii", "mboi")):
        return "DpnII"
    return "DpnII,HinfI"


# ------------------------------------------------------------------ reads entries
def parse_reads_field(value, sep=r"[;,]"):
    """'ERR1;ERR2' or 'r1.fq.gz:r2.fq.gz,single.fq.gz' -> list of read entries (dicts)."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    out = []
    for tok in re.split(sep, str(value)):
        tok = tok.strip()
        if not tok or tok in (".", "nan", "None"):
            continue
        files = [f for f in tok.split(":") if f]
        if all(os.path.exists(f) for f in files):
            stem = re.sub(r"(_R?[12])?(\.f(ast)?q)?(\.gz)?$", "", os.path.basename(files[0]))
            out.append(dict(id=slug(stem), files=[os.path.abspath(f) for f in files],
                            layout="PAIRED" if len(files) == 2 else "SINGLE"))
        elif len(files) == 1:
            out.append(dict(id=files[0], files=[], layout=""))
        else:
            raise ValueError(f"reads '{tok}': local files not found")
    return out


def _split(x):
    if x is None or (isinstance(x, float) and pd.isna(x)):
        return []
    return [t.strip() for t in re.split(r"[;,]", str(x)) if t.strip() not in ("", ".", "nan")]


# ------------------------------------------------------------------ building the table
def _group_rows(rows):
    """Merge rows sharing (name, ref): union of WGS and Hi-C runs (DToL tables have one row
    per WGS x Hi-C combination)."""
    units = {}
    for r in rows:
        key = (r["name"], r["ref"])
        u = units.setdefault(key, dict(name=r["name"], ref=r["ref"], chroms=r.get("chroms", ""),
                                       enzymes=r.get("enzymes", ""), wgs={}, hic={}))
        for kind in ("wgs", "hic"):
            for e in r[kind]:
                u[kind].setdefault(e["id"], e)
        if r.get("chroms"):
            u["chroms"] = r["chroms"]
        if r.get("enzymes"):
            u["enzymes"] = r["enzymes"]
    return list(units.values())


def _attach(entries, layouts, protocols, techs, kind, platforms=()):
    for i, e in enumerate(entries):
        pick = lambda lst: lst[i] if i < len(lst) else (lst[-1] if lst else "")
        e["layout"] = e.get("layout") or pick(layouts)
        e["protocol"] = pick(protocols)
        if kind == "wgs":
            # protocol AND platform: DToL runs with protocol 'WGS' are only typed by their platform
            e["tech"] = pick(techs) or infer_tech(e["protocol"], pick(list(platforms)), e["layout"])
            if e["tech"] in ("", "auto") and e["files"]:
                e["tech"] = guess_tech_from_reads(e["files"][0], e["layout"])
    return entries


def _col(row, *names):
    for n in names:
        if n in row and not (isinstance(row[n], float) and pd.isna(row[n])) and str(row[n]) not in ("", "nan"):
            return str(row[n])
    return ""


def rows_from_table(path):
    sep = "\t" if path.endswith((".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep, dtype=str)
    df = df.loc[:, [c for c in df.columns if not str(c).startswith("Unnamed")]]
    missing = {"name", "ref"} - set(df.columns)
    if missing:
        raise ValueError(f"table {path}: missing mandatory column(s) {sorted(missing)}")
    rows = []
    for _, row in df.iterrows():
        row = row.to_dict()
        # old DToL formatter put study accessions in wgs/hic and runs in run_acession_*
        wgs_field = _col(row, "run_accession_WGS", "run_acession_WGS") if not is_sra(_col(row, "wgs").split(";")[0]) else _col(row, "wgs")
        hic_field = _col(row, "run_accession_Hi-C", "run_acession_Hi-C") if not is_sra(_col(row, "hic").split(";")[0]) else _col(row, "hic")
        wgs_field = wgs_field or _col(row, "wgs")
        hic_field = hic_field or _col(row, "hic")
        wgs = _attach(parse_reads_field(wgs_field), _split(_col(row, "layout_WGS", "wgs_layout")),
                      _split(_col(row, "protocol_WGS", "wgs_protocol")),
                      _split(_col(row, "wgs_tech", "tech")), "wgs", _split(_col(row, "platform_WGS")))
        hic = _attach(parse_reads_field(hic_field), _split(_col(row, "layout_Hi-C", "hic_layout")),
                      _split(_col(row, "protocol_Hi-C", "hic_protocol")), [], "hic")
        for ref in _split(_col(row, "ref")) or [""]:
            rows.append(dict(name=_col(row, "name", "species"), ref=ref, wgs=wgs, hic=hic,
                             chroms=_col(row, "chroms"), enzymes=_col(row, "enzymes")))
    return rows


def rows_from_args(name, ref, wgs, hic, wgs_tech=None, chroms=None, enzymes=None, hic_protocol=None):
    if not name:
        raise ValueError("--name is required when no --table is given")
    wgs_e = _attach(parse_reads_field(wgs), [], [], _split(wgs_tech) if wgs_tech else [], "wgs")
    hic_e = _attach(parse_reads_field(hic), [], _split(hic_protocol) if hic_protocol else [], [], "hic")
    refs = _split(ref) or [""]
    return [dict(name=name, ref=r, wgs=wgs_e, hic=hic_e, chroms=chroms or "", enzymes=enzymes or "") for r in refs]


def to_frame(rows):
    out = []
    for u in _group_rows(rows):
        wgs, hic = list(u["wgs"].values()), list(u["hic"].values())
        hic_protocol = next((e.get("protocol") for e in hic if e.get("protocol")), "")
        species = slug(u["name"])
        rs = ref_slug(u["ref"]) if u["ref"] else "no_ref"
        out.append(dict(
            unit=f"{species}/{rs}", name=u["name"], species=species, ref=u["ref"], ref_slug=rs,
            chroms=os.path.abspath(u["chroms"]) if u["chroms"] else "",
            wgs=";".join(e["id"] for e in wgs),
            wgs_files=";".join(":".join(e["files"]) for e in wgs),
            wgs_layout=";".join(e.get("layout", "") for e in wgs),
            wgs_tech=";".join(e.get("tech", "auto") for e in wgs),
            wgs_protocol=";".join(e.get("protocol", "") for e in wgs),
            hic=";".join(e["id"] for e in hic),
            hic_files=";".join(":".join(e["files"]) for e in hic),
            hic_layout=";".join(e.get("layout", "") for e in hic),
            hic_protocol=hic_protocol,
            enzymes=u["enzymes"] or infer_enzymes(hic_protocol)))
    df = pd.DataFrame(out, columns=SAMPLE_COLUMNS)
    if df.unit.duplicated().any():
        raise ValueError("duplicated (species, reference) units: " + ", ".join(df.unit[df.unit.duplicated()]))
    return df


def guess_tech_from_reads(path, layout=""):
    """Long (>1 kb on average) single-end reads are assumed HiFi, otherwise short reads."""
    import gzip
    try:
        op = gzip.open if path.endswith(".gz") else open
        n, tot = 0, 0
        with op(path, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    n += 1; tot += len(line.strip())
                if n >= 2000:
                    break
        mean = tot / max(n, 1)
        return "hifi" if mean > 1000 else "short"
    except Exception:
        return "short" if layout == "PAIRED" else "auto"


# ------------------------------------------------------------------ access helpers
def reads_of(row, kind):
    """List of read entries (dict id/files/layout/tech/protocol) for a samples.tsv row."""
    ids = _split_keep(row[kind])
    files = _split_keep(row[f"{kind}_files"], n=len(ids))
    layouts = _split_keep(row[f"{kind}_layout"], n=len(ids))
    techs = _split_keep(row["wgs_tech"], n=len(ids)) if kind == "wgs" else [""] * len(ids)
    out = []
    for i, rid in enumerate(ids):
        out.append(dict(id=rid, files=[f for f in files[i].split(":") if f], layout=layouts[i], tech=techs[i] or "auto"))
    return out


def _split_keep(x, n=None):
    if x is None or (isinstance(x, float) and pd.isna(x)) or str(x) == "nan":
        x = ""
    parts = str(x).split(";") if str(x) != "" else []
    if n is not None:
        parts = (parts + [""] * n)[:n]
    return parts


def load_samples(outdir):
    p = os.path.join(outdir, SAMPLES_NAME)
    if not os.path.isfile(p):
        raise FileNotFoundError(f"{p} not found: give --table or --name/--ref/--wgs/--hic")
    return pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)


def resolve_samples(outdir, table=None, name=None, ref=None, wgs=None, hic=None, wgs_tech=None,
                    chroms=None, enzymes=None, hic_protocol=None, units=None):
    """Builds (or reloads) samples.tsv and returns the rows selected by --unit."""
    os.makedirs(outdir, exist_ok=True)
    if table:
        df = to_frame(rows_from_table(table))
    elif name or ref or wgs or hic:
        df = to_frame(rows_from_args(name, ref, wgs, hic, wgs_tech, chroms, enzymes, hic_protocol))
    else:
        df = load_samples(outdir)
    if table or name or ref or wgs or hic:
        path = os.path.join(outdir, SAMPLES_NAME)
        if os.path.isfile(path):  # keep previously declared units, update/add the new ones
            old = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
            df = pd.concat([old[~old.unit.isin(df.unit)], df], ignore_index=True)
        df.to_csv(path, sep="\t", index=False)
        if table or name:
            df = df[df.unit.isin(to_frame(rows_from_table(table)).unit if table else
                                 to_frame(rows_from_args(name, ref, wgs, hic, wgs_tech, chroms, enzymes, hic_protocol)).unit)]
    return select_units(df, units)


def select_units(df, units):
    if not units:
        return df.reset_index(drop=True)
    keep = df.unit.isin(units) | df.species.isin([slug(u) for u in units]) | df.name.isin(units)
    if not keep.any():
        raise ValueError(f"--unit {units}: no matching unit (known: {', '.join(df.unit)})")
    return df[keep].reset_index(drop=True)
