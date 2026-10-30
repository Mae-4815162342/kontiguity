"""Runs 'kontiguity classify' directly on existing .mcool files (no read processing):
builds the minimal unit tree (samples.tsv, genome_chromosomes.tsv, hic/<unit>.mcool) and calls the CLI.

    python tests/real_classify.py OUTDIR NAME file.mcool [NAME file.mcool ...]

Host chromosomes = sequences >= 100 kb that are not mitochondria; every other sequence is a candidate.
"""
import os
import re
import subprocess
import sys

import h5py
import pandas as pd

args = sys.argv[1:sys.argv.index("--")] if "--" in sys.argv else sys.argv[1:]
extra = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []
out = os.path.abspath(args[0])
pairs = list(zip(args[1::2], args[2::2]))
rows = []
for name, mcool in pairs:
    with h5py.File(mcool, "r") as f:
        g = f["resolutions"]["5000"]
        names = [x.decode() for x in g["chroms/name"][:]]
        lens = g["chroms/length"][:]
    unit_dir = os.path.join(out, name, "real")
    os.makedirs(os.path.join(unit_dir, "contigs"), exist_ok=True)
    os.makedirs(os.path.join(unit_dir, "hic"), exist_ok=True)
    mito = re.compile(r"^(Mito|NC_001224\.1|NC_000895\.1)$")
    t = pd.DataFrame(dict(id=names, length=lens))
    t["sequence_type"] = ["organelle" if mito.match(n) else "chromosome" if l >= 100_000 else
                          ("recovered_contig" if re.match(r"^k\d+_", n) else "unplaced") for n, l in zip(names, lens)]
    t["sequence_name"] = t.id
    t["origin"] = "reference"
    t.to_csv(os.path.join(unit_dir, "contigs", "genome_chromosomes.tsv"), sep="\t", index=False)
    dst = os.path.join(unit_dir, "hic", f"{name}.real.mcool")
    if not os.path.exists(dst):
        os.symlink(os.path.abspath(mcool), dst)
    rows.append(dict(unit=f"{name}/real", name=name, species=name, ref="real", ref_slug="real", chroms="", wgs="",
                     wgs_files="", wgs_layout="", wgs_tech="", wgs_protocol="", hic="", hic_files="", hic_layout="",
                     hic_protocol="", enzymes=""))
pd.DataFrame(rows).to_csv(os.path.join(out, "samples.tsv"), sep="\t", index=False)
sys.exit(subprocess.call(["kontiguity", "classify", "-o", out, "--force"] + extra))
