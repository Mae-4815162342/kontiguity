"""Step 'retrieve' (one job per species x reference): candidate contigs.

For every WGS run of the species:
  * Logan contigs exist  -> they are used directly (no reads, no alignment, no assembly);
  * otherwise            -> reads are aligned on the reference, the unmapped reads ("read bin")
                            are pooled per technology and assembled de novo:
                              short reads (Illumina, 10x)  -> MEGAHIT
                              PacBio HiFi / CLR, ONT       -> Flye --meta (reports circular contigs)
Output: retrieve/contigs_raw.fa + contigs_raw.tsv (source, coverage, circularity of every contig)
and retrieve/stats.json (sequencing depth of the host genome, used for relative copy numbers).
"""
import glob
import json
import os
import shutil
import time

import pandas as pd

from kontiguity import logan
from kontiguity.fasta import read_fasta, write_record, gc
from kontiguity.load import read_manifest
from kontiguity.utils.runner import run, require, StepError

PRESET = {"short": "sr", "hifi": "map-hifi", "clr": "map-pb", "ont": "map-ont"}
FLYE_MODE = {"hifi": "--pacbio-hifi", "clr": "--pacbio-raw", "ont": "--nano-hq"}
MIN_READS = {"short": 2000, "hifi": 20, "clr": 50, "ont": 50}


def _tech(m):
    t = m["tech"] if m["tech"] not in ("", "auto") else ("short" if m["layout"] == "PAIRED" else "hifi")
    return t if t in PRESET else "short"


def _wait_stats(path, timeout=600):
    """samtools stats runs in a process substitution: wait until its summary is complete."""
    t0 = time.time()
    while time.time() - t0 < timeout:
        if os.path.isfile(path):
            with open(path) as fh:
                txt = fh.read()
            if "\nRL\t" in txt or "\nCOV\t" in txt or ("SN\tbases mapped (cigar)" in txt and "\nGCF\t" in txt):
                return txt
        time.sleep(2)
    return open(path).read() if os.path.isfile(path) else ""


def _mapped_bases(stats_txt):
    for line in stats_txt.splitlines():
        if line.startswith("SN\tbases mapped (cigar):"):
            return int(line.split("\t")[2])
    return 0


def _count_fastq(path):
    import gzip
    n = 0
    with gzip.open(path, "rt") as fh:
        for n, _ in enumerate(fh, 1):
            pass
    return n // 4


def align_unmapped(lay, m, tech, cfg, logger):
    """Aligns one run on the reference; returns (unmapped fastq files, mapped bases)."""
    log, T = lay.log("retrieve"), cfg["threads"]
    d = os.path.join(lay.retrieve_work, "unmapped")
    os.makedirs(d, exist_ok=True)
    files = [f for f in m["files"].split(":") if f]
    paired = len(files) == 2
    out = [os.path.join(d, f"{m['run']}_1.fq.gz"), os.path.join(d, f"{m['run']}_2.fq.gz")] if paired \
        else [os.path.join(d, f"{m['run']}.fq.gz")]
    stats = os.path.join(d, f"{m['run']}.samtools_stats.txt")
    ok = os.path.join(d, f"{m['run']}.done")   # written only once the read bin of this run is complete
    if os.path.isfile(ok) and all(os.path.isfile(f) for f in out):
        return out, _mapped_bases(open(stats).read())
    aligner = cfg["retrieve"]["aligner"]
    if aligner == "bowtie2" and tech == "short":
        require("bowtie2", "samtools")
        idx = os.path.join(lay.retrieve_work, "bt2", "ref")
        if not os.path.isfile(idx + ".1.bt2") and not os.path.isfile(idx + ".1.bt2l"):
            os.makedirs(os.path.dirname(idx), exist_ok=True)
            run(f"bowtie2-build --threads {T} -q {lay.ref_fa} {idx}", log, logger)
        reads = f"-1 {files[0]} -2 {files[1]} --un-conc-gz {d}/{m['run']}_%.fq.gz" if paired \
            else f"-U {files[0]} --un-gz {out[0]}"
        # paired: bowtie2 writes <run>_1.fq.gz / <run>_2.fq.gz through the % placeholder
        run(f"bowtie2 -p {T} -x {idx} {reads} | samtools stats -@ 2 - > {stats}", log, logger)
        open(ok, "w").close()
        return out, _mapped_bases(open(stats).read())
    require("minimap2", "samtools")
    preset = PRESET[tech]
    os.makedirs(lay.retrieve_work, exist_ok=True)
    mmi = os.path.join(lay.retrieve_work, f"ref.{preset}.mmi")
    if not os.path.isfile(mmi):
        run(f"minimap2 -t {T} -x {preset} -d {mmi}.part {lay.ref_fa} && mv {mmi}.part {mmi}", log, logger)
    if paired:
        fq = f"samtools fastq -f 12 -1 {out[0]} -2 {out[1]} -0 /dev/null -s /dev/null -"
    elif tech == "short":
        fq = f"samtools fastq -f 4 - | gzip -1 > {out[0]}"
    else:
        # long reads: unmapped reads AND reads with < max_aligned of their length aligned
        # (a read spanning a sequence absent from the reference is clipped, not unmapped)
        fq = (f"samtools view - | awk -v F={cfg['retrieve']['long_read_max_aligned']} " + r"""'BEGIN{OFS="\n"}
             { keep = int($2 / 4) % 2; if (!keep) { c=$6; a=0; while (match(c, /[0-9]+[MI=X]/)) { a += substr(c, RSTART, RLENGTH-1); c = substr(c, RSTART+RLENGTH) }
               keep = (a < F * length($10)) }
               if (keep) print "@"$1, $10, "+", $11 }' """ + f"| gzip -1 > {out[0]}")
    run(f"minimap2 -t {T} -ax {preset} {mmi} {' '.join(files)} 2>> {log}.minimap2 "
        f"| samtools view -u -F 0x900 - | tee >(samtools stats - > {stats}) | {fq}", log, logger)
    txt = _wait_stats(stats)
    open(ok, "w").close()
    return out, _mapped_bases(txt)


def assemble(lay, tech, unmapped, cfg, logger):
    """De novo assembly of the pooled unmapped reads of one technology -> (fasta, info dict)."""
    log, T = lay.log("retrieve"), cfg["threads"]
    adir = os.path.join(lay.retrieve_work, f"assembly_{tech}")
    if tech == "short":
        require("megahit")
        final = os.path.join(adir, "final.contigs.fa")
        if not os.path.isfile(final):
            if os.path.isdir(adir):
                shutil.rmtree(adir)
            p1 = [u[0] for u in unmapped if len(u) == 2]
            p2 = [u[1] for u in unmapped if len(u) == 2]
            se = [u[0] for u in unmapped if len(u) == 1]
            reads = (f"-1 {','.join(p1)} -2 {','.join(p2)} " if p1 else "") + (f"-r {','.join(se)}" if se else "")
            run(f"megahit -t {T} {reads} -o {adir} --min-contig-len {cfg['retrieve']['min_contig_len']} "
                f"{cfg['retrieve']['megahit_opts']}", log, logger)
        info = {}
        for h, _ in read_fasta(final):
            f = dict(x.split("=", 1) for x in h.split()[1:] if "=" in x)
            info[h.split()[0]] = dict(coverage=float(f.get("multi", "nan")), coverage_type="megahit_kmer_multiplicity",
                                      circular="unknown", circular_evidence="")
        return final, info, "megahit"
    require("flye")
    final = os.path.join(adir, "assembly.fasta")
    if not os.path.isfile(final):
        reads = " ".join(u[0] for u in unmapped)
        rc = run(f"flye {FLYE_MODE[tech]} {reads} --meta -t {T} -o {adir} {cfg['retrieve']['flye_opts']}",
                 log, logger, check=False)
        if rc != 0 or not os.path.isfile(final):
            logger.warning(f"[retrieve:{lay.unit}] Flye produced no assembly for {tech} reads (too few reads?)")
            return None, {}, "flye"
    info = {}
    ai = os.path.join(adir, "assembly_info.txt")
    if os.path.isfile(ai):
        t = pd.read_csv(ai, sep="\t")
        for r in t.itertuples(index=False):
            info[r[0]] = dict(coverage=float(r[2]), coverage_type="flye_read_depth",
                              circular="yes" if str(r[3]).upper() == "Y" else "no",
                              circular_evidence="Flye circular contig" if str(r[3]).upper() == "Y" else "")
    return final, info, "flye"


def retrieve_unit(lay, cfg, logger):
    os.makedirs(lay.retrieve_dir, exist_ok=True)
    man = read_manifest(lay)
    wgs = man[man.kind == "wgs"]
    if wgs.empty:
        raise StepError(f"{lay.unit}: no WGS run")
    min_len = cfg["retrieve"]["min_contig_len"]
    rows = []
    stats = dict(logan_runs=[], read_runs=[], mapped_bases={}, ref_length=0, assemblies={})
    stats["ref_length"] = sum(len(s) for _, s in read_fasta(lay.ref_fa))
    by_tech = {}
    with open(lay.raw_contigs + ".part", "w") as out:
        for m in wgs.to_dict("records"):
            lf = m["logan"]
            has_logan = bool(lf) and os.path.isfile(lf) and cfg["logan"]["enabled"]
            if has_logan:
                n = logan.extract(lf, out, min_len, f"logan_{m['run']}", rows)
                stats["logan_runs"].append(m["run"])
                logger.info(f"[retrieve:{lay.unit}] {m['run']}: {n} Logan contigs >= {min_len} bp")
            if not m["files"]:
                if has_logan:
                    continue
                if m["source"] == "logan_only":
                    raise StepError(f"{lay.unit}: run {m['run']} was only fetched from Logan and Logan is now disabled: "
                                    f"run 'kontiguity load --no-logan --force' to download its reads")
                raise StepError(f"{lay.unit}: run {m['run']} has neither Logan contigs nor reads")
            tech = _tech(m)
            un, mapped = align_unmapped(lay, m, tech, cfg, logger)
            stats["read_runs"].append(m["run"])
            stats["mapped_bases"][tech] = stats["mapped_bases"].get(tech, 0) + mapped
            by_tech.setdefault(tech, []).append(un)
        for tech, unmapped in by_tech.items():
            if tech == "hifi" and cfg["retrieve"]["read_circles"]:
                from kontiguity.readcircles import find_circles
                circles = find_circles([u[0] for u in unmapped], cfg["retrieve"]["read_circle_min_len"])
                for j, cc in enumerate(circles):
                    cid = f"readcircle_{tech}_{j}"
                    write_record(out, cid, cc["seq"])
                    rows.append(dict(id=cid, source=f"readcircle_{tech}", original_name=f"period_{cc['period']}",
                                     length=len(cc["seq"]), gc=round(gc(cc["seq"]), 4), coverage=cc["copies"],
                                     coverage_type="read_copies", circular="yes",
                                     circular_evidence=f"tandem copies in {cc['support']} long reads (period {cc['period']} bp)"))
                logger.info(f"[retrieve:{lay.unit}] {len(circles)} circular units found as tandem copies in {tech} reads")
            n_reads = sum(_count_fastq(u[0]) for u in unmapped)
            if n_reads < MIN_READS[tech]:
                logger.warning(f"[retrieve:{lay.unit}] only {n_reads} unmapped {tech} reads: no assembly")
                continue
            fa, info, tool = assemble(lay, tech, unmapped, cfg, logger)
            stats["assemblies"][tech] = dict(tool=tool, unmapped_reads=n_reads, fasta=fa)
            if not fa:
                continue
            k = 0
            for h, seq in read_fasta(fa):
                if len(seq) < min_len:
                    continue
                name = h.split()[0]
                cid = f"{tool}_{tech}_{k}"
                k += 1
                write_record(out, cid, seq)
                meta = info.get(name, dict(coverage=float("nan"), coverage_type="", circular="unknown", circular_evidence=""))
                rows.append(dict(id=cid, source=f"{tool}_{tech}", original_name=name, length=len(seq), gc=round(gc(seq), 4),
                                 **meta))
            logger.info(f"[retrieve:{lay.unit}] {tool} ({tech}): {k} contigs >= {min_len} bp from {n_reads} unmapped reads")
    os.replace(lay.raw_contigs + ".part", lay.raw_contigs)
    pd.DataFrame(rows, columns=["id", "source", "original_name", "length", "gc", "coverage", "coverage_type",
                                "circular", "circular_evidence"]).to_csv(lay.raw_info, sep="\t", index=False)
    stats["host_read_depth"] = {t: b / stats["ref_length"] for t, b in stats["mapped_bases"].items() if stats["ref_length"]}
    with open(os.path.join(lay.retrieve_dir, "stats.json"), "w") as fh:
        json.dump(stats, fh, indent=1)
    if not cfg["keep_tmp"]:
        for p in glob.glob(os.path.join(lay.retrieve_work, "unmapped")) + glob.glob(os.path.join(lay.retrieve_work, "*.mmi")) \
                + glob.glob(os.path.join(lay.retrieve_work, "bt2")):
            shutil.rmtree(p, ignore_errors=True) if os.path.isdir(p) else os.remove(p)
        for a in glob.glob(os.path.join(lay.retrieve_work, "assembly_*")):
            for sub in ("intermediate_contigs", "00-assembly", "10-consensus", "20-repeat", "30-contigger", "40-polishing"):
                shutil.rmtree(os.path.join(a, sub), ignore_errors=True)
    logger.info(f"[retrieve:{lay.unit}] {len(rows)} raw contigs")
    return len(rows)
