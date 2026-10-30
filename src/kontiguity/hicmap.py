"""Step 'map' (one job per species x reference): Hi-C contact map on the new reference.

hicstuff aligns the Hi-C reads on contigs/genome.fa (reference + recovered contigs) and produces the
deduplicated pairs; cooler builds the base-resolution matrix and zoomifies it into a balanced .mcool.
(hicstuff's own matrix export is not used: cooler is driven directly for full control of resolutions.)
"""
import glob
import os
import shutil

from kontiguity.load import read_manifest
from kontiguity.utils.runner import run, require, StepError


def _hic_files(lay):
    man = read_manifest(lay)
    hic = man[man.kind == "hic"]
    r1, r2 = [], []
    for files in hic.files:
        f = [x for x in files.split(":") if x]
        if len(f) == 2:
            r1.append(f[0]); r2.append(f[1])
        elif len(f) == 1:
            raise StepError(f"Hi-C run {f[0]} is single-end: hicstuff needs paired reads")
    if not r1:
        raise StepError(f"{lay.unit}: no paired Hi-C reads (run 'kontiguity load')")
    return r1, r2


def map_unit(lay, cfg, logger):
    require("hicstuff", "cooler", "bowtie2")
    m, T, log = cfg["map"], cfg["threads"], lay.log("map")
    if not os.path.isfile(lay.genome_fa):
        raise StepError(f"{lay.genome_fa} missing: run 'kontiguity filter' first")
    os.makedirs(lay.hic_dir, exist_ok=True)
    work = lay.hic_work
    os.makedirs(work, exist_ok=True)
    prefix = f"{lay.row['species']}.{lay.row['ref_slug']}"
    r1, r2 = _hic_files(lay)
    if len(r1) > 1:   # several Hi-C runs: pooled (concatenated gzip streams are valid gzip)
        cat1, cat2 = os.path.join(work, "pooled_R1.fq.gz"), os.path.join(work, "pooled_R2.fq.gz")
        os.makedirs(work, exist_ok=True)
        if not os.path.isfile(cat2):
            run(f"cat {' '.join(r1)} > {cat1} && cat {' '.join(r2)} > {cat2}", log, logger)
        r1, r2 = [cat1], [cat2]
    enz = m["enzymes"] if m["enzymes"] not in ("auto", "", None) else (lay.row.get("enzymes") or "DpnII,HinfI")
    mapping = m["mapping"]
    if mapping == "cutsite" and (enz in ("dnase", "mnase") or str(enz).isdigit()):
        mapping = "normal"   # no ligation motif to split reads on
    pairs = sorted(glob.glob(os.path.join(work, f"{prefix}.valid_idx_pcrfree.pairs*")) +
                   glob.glob(os.path.join(work, f"{prefix}.valid_idx_filtered.pairs*")))
    if not pairs:
        flags = (" -D" if m["duplicates"] else "") + (" -f" if m["filter_events"] else "")
        run(f"hicstuff pipeline -g {lay.genome_fa} -o {work} -P {prefix} -t {T} -a {m['aligner']} -e {enz} -m {mapping} "
            f"-q {m['min_mapq']} -d -n -z False -T {work}/tmp{flags} -F {r1[0]} {r2[0]}", log, logger)
    pairs = [p for p in [os.path.join(work, f"{prefix}.valid_idx_pcrfree.pairs.gz"),
                         os.path.join(work, f"{prefix}.valid_idx_pcrfree.pairs")] +
             sorted(glob.glob(os.path.join(work, "**", f"{prefix}.valid_idx_pcrfree.pairs*"), recursive=True)) +
             sorted(glob.glob(os.path.join(work, "**", f"{prefix}.valid_idx_filtered.pairs*"), recursive=True)) +
             sorted(glob.glob(os.path.join(work, "**", f"{prefix}.valid_idx.pairs*"), recursive=True))
             if os.path.isfile(p)]
    if not pairs:
        raise StepError(f"hicstuff produced no pairs file in {work}")
    pairs = pairs[0]
    cs = os.path.join(lay.hic_dir, "chrom.sizes")
    cat = "zcat" if pairs.endswith(".gz") else "cat"
    run(f"{cat} {pairs} | grep '^#chromsize' | awk '{{print $2\"\\t\"$3}}' > {cs}", log, logger)
    res = sorted(int(x) for x in m["resolutions"])
    base = res[0]
    cool = os.path.join(work, f"{prefix}.{base}.cool")
    run(f"cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 {cs}:{base} {pairs} {cool}", log, logger)
    zres = ",".join(str(r) for r in res if r % base == 0)
    store = lay.mcool_store
    os.makedirs(os.path.dirname(store), exist_ok=True)
    rc = run(f"cooler zoomify -p {T} -r {zres} {'--balance' if m['balance'] else ''} -o {store}.part {cool} "
             f"&& mv {store}.part {store}", log, logger, check=False)
    if rc != 0:   # balancing can fail on very sparse maps: keep raw counts (classification uses raw counts)
        logger.warning(f"[map:{lay.unit}] balanced zoomify failed, writing an unbalanced mcool")
        run(f"cooler zoomify -p {T} -r {zres} -o {store}.part {cool} && mv {store}.part {store}", log, logger)
    if store != lay.mcool:   # --tmp-dir: the output tree links to the map stored in the temporary space
        if os.path.lexists(lay.mcool):
            os.remove(lay.mcool)
        os.symlink(store, lay.mcool)
    for f in glob.glob(os.path.join(work, f"{prefix}*.log")) + glob.glob(os.path.join(work, f"{prefix}.distance_law.txt")):
        shutil.copy(f, lay.hic_dir)
    if not cfg["keep_tmp"]:
        shutil.rmtree(os.path.join(work, "tmp"), ignore_errors=True)
        for f in glob.glob(os.path.join(work, "pooled_R*.fq.gz")):
            os.remove(f)
        os.remove(cool)
    logger.info(f"[map:{lay.unit}] contact map {lay.mcool} (resolutions {zres})")
    return lay.mcool
