"""Small circular molecules from long accurate reads (PacBio HiFi).

A read longer than a circular molecule contains several tandem copies of it: the read is periodic.
For every read of the read bin, the smallest period d >= min_period such that the read matches itself
shifted by d (Hamming identity >= min_identity over the overlap) gives one copy of the circle.
Copies are clustered by period; each cluster with >= min_support reads yields a circular contig.
(Assemblers such as Flye often miss circles shorter than the reads.)
"""
import gzip
from collections import defaultdict


def _read_fastq(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        while True:
            h = fh.readline()
            if not h:
                return
            s = fh.readline().strip()
            fh.readline(); fh.readline()
            yield h[1:].split()[0], s


def period(seq, min_period=1000, seed_len=24, offset=100, check=2000, min_identity=0.98, max_tries=50):
    n = len(seq)
    if n < 2 * min_period + offset:
        return None
    seed = seq[offset:offset + seed_len]
    if "N" in seed:
        return None
    start = offset + min_period
    for _ in range(max_tries):
        pos = seq.find(seed, start)
        if pos < 0:
            return None
        d = pos - offset
        L = min(check, n - offset - d)
        if L >= min(500, d):
            a, b = seq[offset:offset + L], seq[offset + d:offset + d + L]
            mism = sum(1 for x, y in zip(a, b) if x != y)
            if 1 - mism / L >= min_identity:
                return d
        start = pos + 1
    return None


def find_circles(fastqs, min_period=1000, max_period=200_000, min_support=2, tolerance=0.02, max_per_cluster=3,
                 min_identity=0.98):
    """-> list of dict(seq, period, support)."""
    units = []
    for fq in fastqs:
        for name, seq in _read_fastq(fq):
            d = period(seq.upper(), min_period, min_identity=min_identity)
            if d and d <= max_period:
                units.append((d, seq[100:100 + d], len(seq) / d))
    units.sort()
    clusters, cur = [], []
    for u in units:
        if cur and u[0] > cur[0][0] * (1 + tolerance):
            clusters.append(cur); cur = []
        cur.append(u)
    if cur:
        clusters.append(cur)
    out = []
    for cl in clusters:
        if len(cl) < min_support:
            continue
        cl = sorted(cl, key=lambda u: -u[2])          # reads with the most copies first
        for d, s, copies in cl[:max_per_cluster]:
            out.append(dict(seq=s, period=d, support=len(cl), copies=round(sum(c for _, _, c in cl), 1)))
    return out
