"""Logan (assemblies of the whole SRA, https://github.com/IndexThePlanet/Logan).

Contigs are stored per sequencing run at
    https://s3.amazonaws.com/logan-pub/c/<RUN>/<RUN>.contigs.fa.zst
with headers  >RUN_<n> ka:f:<mean k-mer abundance> L:<s>:<target>:<s> ...  (BCALM2-style links).
A contig linked to itself on the same strand (L:+:<n>:+ on contig <n>) closes a cycle in the
assembly graph: it is reported as circular.
"""
import os
import re
import subprocess

import requests

from kontiguity.fasta import read_fasta, write_record, gc
from kontiguity.samples import is_sra
from kontiguity.utils.runner import run, require

_LINK = re.compile(r"^L:([+-]):(\d+):([+-])$")


def url_for(acc, template):
    return template.format(acc=acc)


def available(acc, template, timeout=30):
    """True / False / None (None = could not tell, e.g. no network)."""
    if not is_sra(acc):
        return False
    url = url_for(acc, template)
    try:
        r = requests.head(url, timeout=timeout, allow_redirects=True)
        return r.status_code == 200
    except requests.RequestException:
        pass
    try:
        p = subprocess.run(["wget", "-q", "--spider", "--tries=2", f"--timeout={timeout}", url])
        return p.returncode == 0
    except FileNotFoundError:
        return None


def download(acc, template, outdir, log, logger=None, tries=5, timeout=120):
    """Downloads and decompresses the Logan contigs of one run -> path of the .fa, or None."""
    require("wget", "zstd")
    os.makedirs(outdir, exist_ok=True)
    fa = os.path.join(outdir, f"{acc}.contigs.fa")
    if os.path.isfile(fa) and os.path.getsize(fa) > 0:
        return fa
    zst = fa + ".zst"
    rc = run(f"wget -q -c --tries={tries} --timeout={timeout} -O {zst}.part '{url_for(acc, template)}' && mv {zst}.part {zst}",
             log, logger, check=False)
    if rc != 0 or not os.path.isfile(zst):
        for f in (zst + ".part", zst):
            if os.path.exists(f):
                os.remove(f)
        return None
    run(f"zstd -d -q -f --rm -o {fa}.part {zst} && mv {fa}.part {fa}", log, logger)
    return fa


def parse_header(header):
    """'SRR1_12 ka:f:3.5 L:+:12:+' -> dict(name, index, ka, links, circular)."""
    fields = header.split()
    name = fields[0]
    idx = name.rsplit("_", 1)[-1]
    ka, links = float("nan"), []
    for f in fields[1:]:
        if f.startswith(("ka:f:", "km:f:")):
            try:
                ka = float(f.split(":")[2])
            except ValueError:
                pass
        else:
            m = _LINK.match(f)
            if m:
                links.append(m.groups())
    circular = any(t == idx and s1 == s2 for s1, t, s2 in links)
    return dict(name=name, index=idx, ka=ka, links=len(links), circular=circular)


def extract(fa, out_fa, min_len, source, writer_info):
    """Copies contigs >= min_len from a Logan fasta into out_fa (open handle) and appends
    metadata rows to writer_info (list). Returns the number of contigs kept."""
    n = 0
    for header, seq in read_fasta(fa):
        if len(seq) < min_len:
            continue
        h = parse_header(header)
        cid = f"{source}_{h['index']}"
        write_record(out_fa, cid, seq)
        writer_info.append(dict(id=cid, source=source, original_name=h["name"], length=len(seq), gc=round(gc(seq), 4),
                                coverage=h["ka"], coverage_type="logan_kmer_abundance",
                                circular="yes" if h["circular"] else "no",
                                circular_evidence="Logan graph self-loop" if h["circular"] else ""))
        n += 1
    return n
