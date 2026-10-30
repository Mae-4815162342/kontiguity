"""Simulated dataset for end-to-end tests of kontiguity.

Truth:
  host     : 4 S. cerevisiae chromosomes + mitochondrion (organelle)
  plasmid  : 6.3 kb random circular sequence, 15 copies, contacts enriched at 25 host hotspots
  fragment : 30 kb of chromosome II removed from the reference (chromosomal piece -> must be rejected)
  bacteria : 150 kb random 'contaminant' genome (contacts only with itself)
Outputs (in OUT): reference.fa, chromosomes.csv, truth.fa, wgs_R1/R2.fq.gz, hic_R1/R2.fq.gz, hifi.fq.gz,
                  logan/SRR0000001.contigs.fa (Logan-like contigs of a 'WGS run')
"""
import gzip
import os
import random
import subprocess
import sys

import numpy as np

SRC, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(1)
random.seed(1)


def read_fa(p):
    seqs, h = {}, None
    for line in open(p):
        line = line.strip()
        if line.startswith(">"):
            h = line[1:].split()[0]; seqs[h] = []
        elif h:
            seqs[h].append(line)
    return {k: "".join(v).upper() for k, v in seqs.items()}


def write_fa(path, d):
    with open(path, "w") as fh:
        for k, s in d.items():
            fh.write(f">{k}\n")
            for i in range(0, len(s), 80):
                fh.write(s[i:i + 80] + "\n")


def rand_seq(n, gc=0.4):
    p = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(rng.choice(list("ACGT"), n, p=p))


src = read_fa(SRC)
chroms = {k: src[k] for k in ["NC_001133.9", "NC_001134.8", "NC_001135.5", "NC_001138.5"]}
mito = src["NC_001224.1"]
plasmid = rand_seq(6300, 0.38)
bacteria = rand_seq(150_000, 0.55)
frag_start, frag_len = 300_000, 30_000
chrII = chroms["NC_001134.8"]

truth = dict(chroms); truth["NC_001224.1"] = mito
truth["PLASMID"] = plasmid; truth["BACTERIA"] = bacteria
write_fa(f"{OUT}/truth.fa", truth)
ref = dict(chroms); ref["NC_001134.8"] = chrII[:frag_start] + chrII[frag_start + frag_len:]; ref["NC_001224.1"] = mito
write_fa(f"{OUT}/reference.fa", ref)
with open(f"{OUT}/chromosomes.csv", "w") as fh:
    fh.write("id,sequence_type,sequence_name\n")
    for i, k in enumerate(chroms, 1):
        fh.write(f"{k},chromosome,chr{i}\n")
    fh.write("NC_001224.1,organelle,mitochondrion\n")

COMP = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(COMP)[::-1]


def circ_slice(s, pos, n):
    return (s * (2 + n // len(s)))[pos:pos + n] if pos + n > len(s) else s[pos:pos + n]


def mutate(s, rate=0.002):
    s = list(s)
    for i in np.where(rng.random(len(s)) < rate)[0]:
        s[i] = "ACGT"[(("ACGT".index(s[i]) if s[i] in "ACGT" else 0) + rng.integers(1, 4)) % 4]
    return "".join(s)


# ---------------------------------------------------------------- WGS (paired 150 bp, ~15x host)
def wgs_pairs(seq, n, circular=False, L=150, ins=400):
    out = []
    for _ in range(n):
        p = int(rng.integers(0, len(seq) - (0 if circular else ins)))
        frag = circ_slice(seq, p, ins) if circular else seq[p:p + ins]
        if len(frag) < ins:
            continue
        r1, r2 = frag[:L], rc(frag[-L:])
        if rng.random() < 0.5:
            r1, r2 = r2, r1
        out.append((mutate(r1), mutate(r2)))
    return out


pairs = []
host_len = sum(len(s) for s in chroms.values())
for s in chroms.values():
    pairs += wgs_pairs(s, int(15 * len(s) / 300))
pairs += wgs_pairs(mito, int(60 * len(mito) / 300))
pairs += wgs_pairs(plasmid, int(15 * 15 * len(plasmid) / 300), circular=True)
pairs += wgs_pairs(bacteria, int(8 * len(bacteria) / 300))
random.shuffle(pairs)
q = "I" * 150
with gzip.open(f"{OUT}/wgs_R1.fq.gz", "wt") as f1, gzip.open(f"{OUT}/wgs_R2.fq.gz", "wt") as f2:
    for i, (a, b) in enumerate(pairs):
        f1.write(f"@w{i}\n{a}\n+\n{q[:len(a)]}\n"); f2.write(f"@w{i}\n{b}\n+\n{q[:len(b)]}\n")
print("WGS pairs", len(pairs))

# ---------------------------------------------------------------- HiFi-like long reads (8 kb, 0.1% err)
with gzip.open(f"{OUT}/hifi.fq.gz", "wt") as fh:
    k = 0
    for name, s, depth, circ in [*[(n, s, 10, False) for n, s in chroms.items()], ("m", mito, 30, False),
                                  ("p", plasmid, 150, True), ("b", bacteria, 10, False)]:
        n = int(depth * len(s) / 8000)
        for _ in range(n):
            L = int(rng.normal(15000, 3000)) if circ else int(min(rng.normal(8000, 1500), len(s) - 1))
            L = max(L, 1500)
            p = int(rng.integers(0, len(s) - (0 if circ else L)))
            r = circ_slice(s, p, L) if circ else s[p:p + L]
            r = mutate(r, 0.001)
            if rng.random() < 0.5:
                r = rc(r)
            fh.write(f"@h{k}\n{r}\n+\n{'I' * len(r)}\n"); k += 1
print("HiFi reads", k)

# ---------------------------------------------------------------- Hi-C (DpnII-like, 100 bp ends)
names = list(chroms)
lens = np.array([len(chroms[n]) for n in names], float)
hot = [(names[int(i)], int(rng.integers(0, int(lens[int(i)]) - 20000))) for i in rng.choice(len(names), 25, p=lens / lens.sum())]


def pos_cis(ch, x):
    L = len(truth[ch])
    d = int(np.exp(rng.uniform(np.log(1000), np.log(L))))
    y = x + d if rng.random() < 0.5 else x - d
    return min(max(y, 0), L - 101)


def rand_host():
    i = int(rng.choice(len(names), p=lens / lens.sum()))
    return names[i], int(rng.integers(0, int(lens[i]) - 101))


def end(ch, p, circular=False):
    s = truth[ch]
    r = circ_slice(s, p, 100) if circular else s[p:p + 100]
    return mutate(r if rng.random() < 0.5 else rc(r))


N = 350_000
with gzip.open(f"{OUT}/hic_R1.fq.gz", "wt") as f1, gzip.open(f"{OUT}/hic_R2.fq.gz", "wt") as f2:
    q = "I" * 100
    for i in range(N):
        u = rng.random()
        if u < 0.035:      # plasmid
            a = ("PLASMID", int(rng.integers(0, len(plasmid))), True)
            v = rng.random()
            if v < 0.15:
                b = ("PLASMID", int(rng.integers(0, len(plasmid))), True)
            elif v < 0.75:
                ch, x = hot[int(rng.integers(0, len(hot)))]
                b = (ch, x + int(rng.integers(0, 19000)), False)
            else:
                b = (*rand_host(), False)
        elif u < 0.06:     # bacteria (own contacts only)
            x = int(rng.integers(0, len(bacteria) - 101))
            y = min(max(x + int(rng.normal(0, 20000)), 0), len(bacteria) - 101)
            a, b = ("BACTERIA", x, False), ("BACTERIA", y, False)
        elif u < 0.07:     # mito: random ligation
            a = ("NC_001224.1", int(rng.integers(0, len(mito) - 101)), False)
            b = (*rand_host(), False) if rng.random() < 0.3 else ("NC_001224.1", int(rng.integers(0, len(mito) - 101)), False)
        else:              # host
            ch, x = rand_host()
            a = (ch, x, False)
            b = (ch, pos_cis(ch, x), False) if rng.random() < 0.8 else (*rand_host(), False)
        r1, r2 = end(*a), end(*b)
        f1.write(f"@c{i}\n{r1}\n+\n{q}\n"); f2.write(f"@c{i}\n{r2}\n+\n{q}\n")
print("Hi-C pairs", N)

# ---------------------------------------------------------------- Logan-like contigs for a fake run
os.makedirs(f"{OUT}/logan/c/SRR0000001", exist_ok=True)
with open(f"{OUT}/logan/SRR0000001.contigs.fa", "w") as fh:
    k = 0
    # chromosomal contigs (abundance ~15), fragment, plasmid (self-loop, ~225), bacteria pieces (~8), mito (~60)
    for name, s in chroms.items():
        for p in range(0, len(s) - 20000, 20000):
            fh.write(f">SRR0000001_{k} ka:f:{rng.normal(15, 2):.1f} L:+:{k + 1}:+\n{s[p:p + 20000]}\n"); k += 1
    fh.write(f">SRR0000001_{k} ka:f:224.7 L:+:{k}:+ L:-:{k}:-\n{plasmid}\n"); k += 1
    for p in range(0, len(bacteria), 30000):
        fh.write(f">SRR0000001_{k} ka:f:{rng.normal(8, 1):.1f}\n{bacteria[p:p + 30000]}\n"); k += 1
    fh.write(f">SRR0000001_{k} ka:f:61.0 L:+:{k}:+\n{mito}\n"); k += 1
    fh.write(f">SRR0000001_{k} ka:f:14.2\n{chrII[frag_start:frag_start + frag_len]}\n"); k += 1
    fh.write(f">SRR0000001_{k} ka:f:3.0\n{rand_seq(600)}\n"); k += 1   # too short
subprocess.run(["zstd", "-q", "-f", f"{OUT}/logan/SRR0000001.contigs.fa", "-o",
                f"{OUT}/logan/c/SRR0000001/SRR0000001.contigs.fa.zst"], check=True)
print("hotspots", hot[:3], "...")
