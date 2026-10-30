"""Small, dependency-free FASTA helpers (gzip aware)."""
import gzip
import re


def _open(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def read_fasta(path):
    """Yields (header_without_>, sequence)."""
    header, chunks = None, []
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line[0] == ">":
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:], []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def write_record(fh, name, seq, width=80):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


def gc(seq):
    s = seq.upper()
    acgt = sum(s.count(b) for b in "ACGT")
    return (s.count("G") + s.count("C")) / acgt if acgt else float("nan")


def lengths(path):
    return {h.split()[0]: len(s) for h, s in read_fasta(path)}


# ---------------------------------------------------------------- ENA headers
# >ENA|OX637620|OX637620.1 Tausonia pullulans genome assembly, chromosome: 1
_ENA_TYPES = [
    ("mitochondri", "organelle"), ("chloroplast", "organelle"), ("plastid", "organelle"),
    ("apicoplast", "organelle"), ("kinetoplast", "organelle"),
    ("unlocalized", "unlocalized"), ("unplaced", "unplaced"),
    ("chromosome", "chromosome"), ("linkage group", "chromosome"),
    ("scaffold", "scaffold"), ("contig", "contig"), ("plasmid", "plasmid"),
]


def parse_ena_header(header):
    """-> (clean_id, sequence_type, sequence_name, description)."""
    first, _, desc = header.partition(" ")
    parts = first.split("|")
    clean = parts[2] if len(parts) >= 3 and parts[0] == "ENA" else first
    d = desc.lower()
    stype = "unknown"
    # the part after the last comma usually reads 'chromosome: 3' / 'mitochondrion'
    tail = d.split(",")[-1].strip()
    for key, t in _ENA_TYPES:
        if key in tail:
            stype = t
            break
    else:
        for key, t in _ENA_TYPES:
            if key in d:
                stype = t
                break
    m = re.search(r"(chromosome|linkage group)\s*:?\s*([\w.-]+)", desc.split(",")[-1].strip(), re.I)
    if m:
        name = m.group(2)
    elif stype == "organelle":
        name = next((o for o in ("mitochondrion", "chloroplast", "plastid", "apicoplast", "kinetoplast") if o[:7] in d), tail)
    else:
        name = clean
    return clean, stype, name, desc
