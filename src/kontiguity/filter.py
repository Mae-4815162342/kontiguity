"""Step 'filter' (one job per species x reference): sorting the candidate contigs and building
the new reference.

1. length filter;
2. re-alignment of every contig on the reference (minimap2 -x asm20): a contig is kept when less than
   `max_ref_cov` of its length aligns with identity >= `min_identity` (i.e. it does not match the genome);
3. redundancy removal between kept contigs (Logan runs of the same species overlap a lot);
4. circularity of contigs not already characterised by Logan/Flye (end overlap, from kontiguity.utils.circularity);
5. relative copy number = contig coverage / host coverage of the same source;
6. new reference: contigs/genome.fa = reference + kept contigs (+ genome_chromosomes.tsv).
"""
import json
import os
from collections import defaultdict

import numpy as np
import pandas as pd

from kontiguity.fasta import read_fasta, write_record
from kontiguity.utils.runner import run, require, StepError

PAF_COLS = ["q", "qlen", "qs", "qe", "strand", "t", "tlen", "ts", "te", "match", "block", "mapq"]


def read_paf(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 12:
                rows.append(f[:12])
    df = pd.DataFrame(rows, columns=PAF_COLS)
    for c in ["qlen", "qs", "qe", "tlen", "ts", "te", "match", "block", "mapq"]:
        df[c] = df[c].astype(int)
    df["identity"] = df["match"] / df["block"].clip(lower=1)
    return df


def covered_fraction(intervals, length):
    if not intervals or length <= 0:
        return 0.0
    iv = sorted(intervals)
    tot, (cs, ce) = 0, iv[0]
    for s, e in iv[1:]:
        if s > ce:
            tot += ce - cs
            cs, ce = s, e
        else:
            ce = max(ce, e)
    tot += ce - cs
    return min(1.0, tot / length)


def ref_coverage(paf, min_identity):
    """Per contig: fraction of its length covered by reference alignments, best hit, best identity."""
    out = {}
    for q, g in paf.groupby("q"):
        ok = g[g.identity >= min_identity]
        frac = covered_fraction(list(zip(ok.qs, ok.qe)), int(g.qlen.iloc[0])) if len(ok) else 0.0
        best = g.loc[g.match.idxmax()]
        out[q] = (frac, best.t, round(float(best.identity), 3))
    return out


def self_circular(seqs, cfg_f, threads):
    """End-overlap circularity (kontiguity.utils.circularity), for contigs of unknown topology."""
    if not seqs:
        return {}
    from kontiguity.utils.circularity import check_self_circularity
    hits = check_self_circularity(seqs, 13, cfg_f["circularity_min_overlap"], cfg_f["circularity_max_overlap"],
                                  cfg_f["circularity_max_mismatch"], processes=threads)
    return {h["sequences"][0]: (h["overlap_len"], f"self-overlap of {h['overlap_len']} bp ({100 * h['mismatch_rate']:.1f}% "
                                                  f"mismatch, trimmed)") for h in hits}


def filter_unit(lay, cfg, logger):
    require("minimap2")
    f, T, log = cfg["filter"], cfg["threads"], lay.log("filter")
    os.makedirs(lay.contigs_dir, exist_ok=True)
    if not os.path.isfile(lay.raw_contigs):
        raise StepError(f"{lay.raw_contigs} missing: run 'kontiguity retrieve' first")
    info = pd.read_csv(lay.raw_info, sep="\t", keep_default_na=False, na_values=[""]).set_index("id")
    tmp = lay.filter_work
    os.makedirs(tmp, exist_ok=True)

    # 1. length
    seqs = {h.split()[0]: s for h, s in read_fasta(lay.raw_contigs) if f["min_contig_len"] <= len(s) <= f["max_contig_len"]}
    sel = os.path.join(tmp, "len_ok.fa")
    with open(sel, "w") as out:
        for k, s in seqs.items():
            write_record(out, k, s)
    logger.info(f"[filter:{lay.unit}] {len(seqs)} contigs within length bounds")

    # 2. re-alignment on the reference
    paf_path = os.path.join(lay.contigs_dir, "contigs_vs_reference.paf")
    if seqs:
        run(f"minimap2 -t {T} -x {f['minimap2_preset']} --secondary=no {lay.ref_fa} {sel} > {paf_path}", log, logger)
        paf = read_paf(paf_path)
    else:
        open(paf_path, "w").close()
        paf = pd.DataFrame(columns=PAF_COLS + ["identity"])
    cov = ref_coverage(paf, f["min_identity"])
    refcov = {k: cov.get(k, (0.0, "", np.nan)) for k in seqs}
    kept = [k for k in seqs if refcov[k][0] < f["max_ref_cov"]]
    logger.info(f"[filter:{lay.unit}] {len(kept)}/{len(seqs)} contigs do not match the reference "
                f"(< {f['max_ref_cov']:.0%} covered at >= {f['min_identity']:.0%} identity)")

    # host coverage per source (Logan: length-weighted median k-mer abundance of contigs matching the genome)
    host_cov = {}
    stats_path = os.path.join(lay.retrieve_dir, "stats.json")
    stats = json.load(open(stats_path)) if os.path.isfile(stats_path) else {}
    for src, g in info.loc[[k for k in seqs if refcov[k][0] >= 0.9]].groupby("source") if seqs else []:
        if src.startswith("logan_") and g.coverage.notna().any():
            v, w = g.coverage.values.astype(float), g.length.values.astype(float)
            o = np.argsort(v)
            cw = np.cumsum(w[o])
            host_cov[src] = float(v[o][np.searchsorted(cw, cw[-1] / 2)])
    for tech, d in (stats.get("host_read_depth") or {}).items():
        for src in info.source.unique():
            if src.endswith("_" + tech) and not src.startswith("logan_"):
                host_cov[src] = d

    # 3. redundancy removal
    dup_of = {}
    if f["dedup"] and len(kept) > 1:
        kfa = os.path.join(tmp, "kept.fa")
        with open(kfa, "w") as out:
            for k in kept:
                write_record(out, k, seqs[k])
        spaf = os.path.join(tmp, "self.paf")
        run(f"minimap2 -t {T} -x {f['minimap2_preset']} -X {kfa} {kfa} > {spaf}", log, logger)
        sp = read_paf(spaf)
        sp = sp[(sp.q != sp.t) & (sp.identity >= f["dedup_identity"])]
        contain = defaultdict(dict)
        for (q, t), g in sp.groupby(["q", "t"]):   # -X reports each pair once: use both sides
            cq = covered_fraction(list(zip(g.qs, g.qe)), int(g.qlen.iloc[0]))
            ct = covered_fraction(list(zip(g.ts, g.te)), int(g.tlen.iloc[0]))
            contain[q][t] = max(contain[q].get(t, 0), cq)
            contain[t][q] = max(contain[t].get(q, 0), ct)
        contain = {k: list(v.items()) for k, v in contain.items()}
        # priority: circular first, assembled consensus before single-read circles, then longest
        prio = sorted(kept, key=lambda k: (info.at[k, "circular"] != "yes", str(info.at[k, "source"]).startswith("readcircle"),
                                           -len(seqs[k]), k))
        keep_set = []
        ks = set()
        for k in prio:
            hit = next((t for t, c in contain.get(k, []) if c >= f["dedup_cov"] and t in ks and len(seqs[t]) >= len(seqs[k]) * 0.98), None)
            if hit is None:
                hit = next((t for t, c in contain.get(k, []) if c >= f["dedup_cov"] and t in ks), None)
            if hit:
                dup_of[k] = hit
            else:
                keep_set.append(k)
                ks.add(k)
        kept = keep_set
        logger.info(f"[filter:{lay.unit}] {len(dup_of)} redundant contigs removed, {len(kept)} left")

    # 4. circularity for contigs of unknown topology
    unknown = {k: seqs[k] for k in kept if str(info.at[k, "circular"]) in ("unknown", "nan", "")}
    circ = self_circular(unknown, f, T)
    for k, (ovl, _) in circ.items():   # remove the duplicated end so that Hi-C reads map uniquely
        seqs[k] = seqs[k][:-ovl]

    # 5. final table and new reference
    n_dups = pd.Series(dup_of).value_counts().to_dict() if dup_of else {}
    rows = []
    with open(lay.contigs_fa + ".part", "w") as out:
        for i, k in enumerate(sorted(kept, key=lambda k: -len(seqs[k])), 1):
            nid = f"ktg_{i:05d}"
            write_record(out, nid, seqs[k])
            r = info.loc[k]
            src = r["source"]
            circular, evidence = r["circular"], r["circular_evidence"]
            if k in circ:
                circular, evidence = "yes", circ[k][1]
            elif k in unknown:
                circular = "no"
            hc = host_cov.get(src, np.nan)
            rows.append(dict(id=nid, raw_id=k, source=src, original_name=r["original_name"], length=len(seqs[k]),
                             gc=r["gc"], coverage=r["coverage"], coverage_type=r["coverage_type"],
                             host_coverage=round(hc, 3) if hc == hc else np.nan,
                             relative_copy_number=round(float(r["coverage"]) / hc, 3) if hc and hc == hc and r["coverage"] == r["coverage"] else np.nan,
                             circular=circular, circular_evidence=evidence,
                             ref_cov_fraction=round(refcov[k][0], 3), ref_best_hit=refcov[k][1], ref_best_identity=refcov[k][2],
                             n_redundant=n_dups.get(k, 0)))
    os.replace(lay.contigs_fa + ".part", lay.contigs_fa)
    df = pd.DataFrame(rows)
    df.to_csv(lay.contigs_info, sep="\t", index=False)
    chroms = pd.read_csv(lay.ref_chroms, sep="\t", keep_default_na=False)
    run(f"cat {lay.ref_fa} {lay.contigs_fa} > {lay.genome_fa}.part && mv {lay.genome_fa}.part {lay.genome_fa}", log, logger)
    extra = pd.DataFrame(dict(id=df.get("id", []), sequence_type="recovered_contig", sequence_name=df.get("raw_id", []),
                              length=df.get("length", []), origin=df.get("source", []), description=""))
    pd.concat([chroms, extra], ignore_index=True).to_csv(lay.genome_chroms, sep="\t", index=False)
    if not cfg["keep_tmp"]:
        import shutil
        shutil.rmtree(tmp, ignore_errors=True)
    logger.info(f"[filter:{lay.unit}] new reference {lay.genome_fa}: reference + {len(df)} contigs "
                f"({int((df.circular == 'yes').sum()) if len(df) else 0} circular)")
    return len(df)
