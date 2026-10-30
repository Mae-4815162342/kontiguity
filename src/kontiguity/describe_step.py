"""Step 'describe' (one job per species x reference): descriptive outputs of the recovered contigs.

* describe/recovered_contigs.tsv  every recovered contig: origin, length, GC, circularity, coverage,
                                  Hi-C contacts and classification call (from the report);
* describe/known_hits.tsv         (--known FASTA) which recovered contigs / reference sequences match
                                  known sequences (e.g. the 2-micron plasmid), with coverage, identity
                                  and the classification of the matching contigs;
* describe/plots/                 the kontiguity descriptive figures (contact maps, signals, mini
                                  matrices; scripts/describe/describe.py) for the selected contigs.
"""
import os

import pandas as pd

from kontiguity.filter import read_paf, covered_fraction
from kontiguity.utils.runner import run, require, StepError


def _recovered_table(lay):
    info = pd.read_csv(lay.contigs_info, sep="\t", keep_default_na=False) if os.path.isfile(lay.contigs_info) else pd.DataFrame()
    if os.path.isfile(lay.report):
        rep = pd.read_csv(lay.report, sep="\t", keep_default_na=False)
        cols = ["id", "call", "p_hitchhiker", "hic_trans_contacts", "molecule", "intra_fraction", "effect_size",
                "best_host_window", "flags"]
        rep = rep[[c for c in cols if c in rep.columns]]
        info = info.merge(rep, on="id", how="left") if len(info) else rep
    return info


def known_hits(lay, known, cfg, logger, out):
    """minimap2 of known sequences against the new reference (reference + recovered contigs)."""
    require("minimap2")
    paf = os.path.join(os.path.dirname(out), "known_vs_genome.paf")
    run(f"minimap2 -t {cfg['threads']} -x asm20 --secondary=yes -N 50 {lay.genome_fa} {known} > {paf}",
        lay.log("describe"), logger)
    p = read_paf(paf)
    rows = []
    if len(p):
        types = pd.read_csv(lay.genome_chroms, sep="\t", keep_default_na=False).set_index("id").sequence_type.to_dict()
        rep = pd.read_csv(lay.report, sep="\t", keep_default_na=False).set_index("id") if os.path.isfile(lay.report) else None
        for (q, t), g in p.groupby(["q", "t"]):
            ok = g[g.identity >= 0.8]
            if ok.empty:
                continue
            rows.append(dict(
                known=q, known_length=int(g.qlen.iloc[0]), hit=t, hit_type=types.get(t, ""), hit_length=int(g.tlen.iloc[0]),
                known_covered=round(covered_fraction(list(zip(ok.qs, ok.qe)), int(g.qlen.iloc[0])), 3),
                hit_covered=round(covered_fraction(list(zip(ok.ts, ok.te)), int(g.tlen.iloc[0])), 3),
                identity=round(float((ok.match.sum()) / max(ok.block.sum(), 1)), 4),
                call=rep.at[t, "call"] if rep is not None and t in rep.index else "",
                p_hitchhiker=rep.at[t, "p_hitchhiker"] if rep is not None and t in rep.index else ""))
    df = pd.DataFrame(rows, columns=["known", "known_length", "hit", "hit_type", "hit_length", "known_covered",
                                     "hit_covered", "identity", "call", "p_hitchhiker"])
    df = df.sort_values(["known", "known_covered"], ascending=[True, False])
    df.to_csv(out, sep="\t", index=False)
    # summary per known sequence: is it in the recovered contigs?
    for k in sorted(set(p.q)) if len(p) else []:
        d = df[(df.known == k) & (df.hit_type == "recovered_contig")]
        logger.info(f"[describe:{lay.unit}] known '{k}': " + (
            "NOT found among the recovered contigs" if d.empty else
            f"found in {len(d)} recovered contig(s) ({', '.join(d.hit)}), up to {d.known_covered.max():.0%} of its length "
            f"in one contig"))
    return df


def describe_unit(lay, cfg, logger):
    c = cfg["describe"]
    d = os.path.join(lay.unit_dir, "describe")
    os.makedirs(d, exist_ok=True)
    if not os.path.isfile(lay.genome_chroms):
        raise StepError(f"{lay.genome_chroms} missing: run 'kontiguity filter' (or map --genome) first")

    rec = _recovered_table(lay)
    rec.to_csv(os.path.join(d, "recovered_contigs.tsv"), sep="\t", index=False)
    logger.info(f"[describe:{lay.unit}] {len(rec)} recovered contigs listed in describe/recovered_contigs.tsv")

    if c.get("known"):
        if not os.path.isfile(lay.genome_fa):
            raise StepError(f"--known needs {lay.genome_fa}")
        known_hits(lay, c["known"], cfg, logger, os.path.join(d, "known_hits.tsv"))

    if not c.get("plots", True) or not os.path.isfile(lay.mcool):
        if not os.path.isfile(lay.mcool):
            logger.warning(f"[describe:{lay.unit}] no contact map: descriptive figures skipped")
        return len(rec)

    # selection of the contigs to draw (each figure costs a few seconds)
    sel = []
    types = pd.read_csv(lay.genome_chroms, sep="\t", keep_default_na=False)
    if c["contigs"] == "all":
        sel = list(types[types.sequence_type != "chromosome"].id)
    elif len(rec):
        r = rec.copy()
        r["_hic"] = pd.to_numeric(r.get("hic_trans_contacts", 0), errors="coerce").fillna(0)
        r = r[r._hic > 0]
        if c["contigs"] == "candidates" and "call" in r:
            r = r[r.call.isin(["STRONG_CANDIDATE", "WEAK_SIGNAL", "QC_FAIL"]) | (r.get("circular", "") == "yes")]
        sel = list(r.sort_values("_hic", ascending=False).id)
    if c.get("known") and os.path.isfile(os.path.join(d, "known_hits.tsv")):
        kh = pd.read_csv(os.path.join(d, "known_hits.tsv"), sep="\t")
        sel = list(dict.fromkeys(list(kh[kh.hit_type == "recovered_contig"].hit) + sel))   # known matches first
    sel = sel[: int(c["max_contigs"])]
    if not sel:
        logger.info(f"[describe:{lay.unit}] no contig with Hi-C contacts to draw")
        return len(rec)

    # chromosome table in the CSV format expected by describe.py (id, sequence_type, sequence_name)
    chroms_csv = os.path.join(d, "chromosomes.csv")
    types[["id", "sequence_type", "sequence_name"]].to_csv(chroms_csv, index=False)
    script = os.path.join(os.path.dirname(__file__), "scripts", "describe", "describe.py")
    index = lay.genome_fa[: -len(".fa")]          # describe.py looks for <index>.fa
    binnings = ";".join(str(b) for b in c["binnings"])
    formats = ",".join(c["formats"])
    plots = os.path.join(d, "plots")
    run(f"python3 {script} '{lay.row['name']}' {lay.unit_dir} {plots} {chroms_csv} . . {cfg['classify']['host_min_len']} "
        f"{','.join(sel)} {index} '{binnings}' . {lay.mcool} {formats} false {'false' if c.get('mini', True) else 'true'} .",
        lay.log("describe"), logger)
    logger.info(f"[describe:{lay.unit}] figures of {len(sel)} contig(s) in describe/plots/")
    return len(rec)
