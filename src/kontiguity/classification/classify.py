"""Step 'classify' (one job per species x reference).

Decision unit = 'molecule':
  * Hi-C communities (Louvain on the sequence x sequence contact graph) group the fragments of one
    molecule; a small community (< community_small bp) without any host chromosome is merged and scored
    as a whole (e.g. a plasmid split into 3 contigs);
  * a contig inside a Mb-scale community that shares >= community_frac of its contacts with other
    contigs of that community is a piece of a genome: 'chromosomal fragment' if the community holds
    host chromosomes, 'other genome' otherwise (endosymbiont, contaminant);
  * the classifier (12 multi-scale features -> MLP ensemble) gives p = P(2-micron / Ddp5-like
    'hitchhiker' contact profile);
  * calls:  STRONG_CANDIDATE  p > threshold, amplitude >= min_effect, intra-contact share in
                              [min_intra, max_intra]
            WEAK_SIGNAL       p > threshold but amplitude < min_effect
            QC_FAIL           p > threshold but no intra contacts / organelle-like
            NOT_EVALUABLE     < min_depth trans contacts
            OTHER             p <= threshold
"""
import json
import os

import networkx as nx
import numpy as np
import pandas as pd

from kontiguity.classification.features import hand_feats, effect_size, autocorr1, smooth_oe, NAMES
from kontiguity.classification.model import MLPEnsemble
from kontiguity.classification.profiles import Profiles
from kontiguity.utils.runner import StepError


def communities(P, seed=0, min_contacts=5):
    C = P.seq_contacts.tocoo()
    T = np.asarray(P.seq_contacts.sum(1)).ravel()
    keep = (C.row < C.col) & (C.data >= min_contacts)
    G = nx.Graph()
    G.add_nodes_from(range(len(P.names)))
    w = C.data[keep] / np.sqrt(T[C.row[keep]] * T[C.col[keep]])
    G.add_weighted_edges_from(zip(C.row[keep], C.col[keep], w))
    lab = np.empty(len(P.names), int)
    for k, cm in enumerate(nx.community.louvain_communities(G, weight="weight", seed=seed)):
        lab[list(cm)] = k
    return lab


def plot_profile(r, P, title, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rng = np.random.default_rng(0)
    y = smooth_oe(r, P.w, P.hmask)
    nulls = np.stack([smooth_oe(rng.multinomial(int(r.sum()), P.w).astype(float), P.w, P.hmask) for _ in range(20)])
    lo, hi = np.percentile(nulls, [2.5, 97.5], 0)
    fig, ax = plt.subplots(figsize=(13, 2.6))
    x = np.arange(len(y))
    ax.fill_between(x, lo, hi, color="#c9c8c0", lw=0, label="null 95% (same depth)")
    ax.plot(x, y, lw=.7, color="#2a78d6")
    for b in np.where(np.diff(P.chrom))[0]:
        ax.axvline(b, color="#999", lw=.3)
    ax.set_xlim(0, len(y)); ax.set_ylim(-3, 3.5)
    ax.set_ylabel("log2 obs/exp"); ax.set_xlabel("host genome (bins, chromosomes concatenated; 50-kb smoothing)")
    ax.set_title(title, loc="left", fontsize=9); ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout(); fig.savefig(path, dpi=110); plt.close(fig)


def classify_unit(lay, cfg, logger):
    c = cfg["classify"]
    if not os.path.isfile(lay.mcool):
        raise StepError(f"{lay.mcool} missing: run 'kontiguity map' first")
    os.makedirs(lay.class_dir, exist_ok=True)
    chroms = pd.read_csv(lay.genome_chroms, sep="\t", keep_default_na=False)
    info = pd.read_csv(lay.contigs_info, sep="\t", keep_default_na=False).set_index("id") if os.path.isfile(lay.contigs_info) else pd.DataFrame()
    host = chroms[chroms.sequence_type.isin(c["host_types"]) & (chroms.length.astype(int) >= c["host_min_len"])].id
    if host.empty:   # untyped reference: long sequences are the host
        host = chroms[(chroms.sequence_type != "recovered_contig") & (chroms.length.astype(int) >= c["host_min_len"])].id
    model = MLPEnsemble(c["model"])
    res = c["resolution"] if c["resolution"] else model.resolution
    P = Profiles(lay.mcool, set(host), res)
    logger.info(f"[classify:{lay.unit}] {P.host.sum()} host chromosomes, {len(P.cand_names)} candidate sequences")

    ctype = chroms.set_index("id").sequence_type.to_dict()
    corigin = chroms.set_index("id").origin.to_dict() if "origin" in chroms else {}
    lab = communities(P, seed=c["seed"])
    seq_i = {n: i for i, n in enumerate(P.names)}
    comm_len = pd.Series(P.lengths).groupby(lab).sum()
    comm_host = pd.Series(P.host).groupby(lab).sum()
    CCm = P.seq_contacts
    members = pd.Series(np.arange(len(P.names))).groupby(lab).apply(list).to_dict()   # community -> sequence indices

    # molecules = merged small communities without host chromosome
    cand = [n for n in P.cand_names if c["include_reference_scaffolds"] or ctype.get(n) == "recovered_contig"]
    molecule = {}
    for n in cand:
        k = lab[seq_i[n]]
        small = comm_len[k] < c["community_small"] and comm_host[k] == 0
        molecule[n] = tuple(sorted(P.names[j] for j in members[k] if P.names[j] in P.index)) if small else (n,)

    cache, rows = {}, []
    thr = c["threshold"]
    for n in cand:
        i = seq_i[n]
        k = lab[i]
        big = not (comm_len[k] < c["community_small"] and comm_host[k] == 0)
        tot_i = P.trans[P.index[n]]
        selfc = P.self_contacts[P.index[n]]
        others = [j for j in members[k] if j != i and not P.host[j]] if big else []
        comm_frac = CCm[i, others].sum() / (CCm[i].sum() - CCm[i, i] + 1) if big and others else 0.0
        mol = molecule[n]
        if mol not in cache:
            r = P.counts(list(mol))
            depth = r.sum()
            p = float(model.predict_proba(hand_feats(r, P.w, P.hmask, P.chrom))[0]) if depth >= c["min_depth"] else np.nan
            if len(mol) > 1:
                ii = [seq_i[m] for m in mol]
                sc = P.self_contacts[[P.index[m] for m in mol]].sum() + CCm[ii][:, ii].sum() / 2 - CCm[ii, ii].sum() / 2
            else:
                sc = selfc
            intra = sc / (sc + depth + 1)
            eff = ac = np.nan
            if p == p and p > thr:
                eff, ac = effect_size(r, P.w, P.hmask, c["n_null"], c["seed"]), autocorr1(r, P.w, P.hmask)
            feats = hand_feats(r, P.w, P.hmask, P.chrom)[0] if depth > 0 else np.full(len(NAMES), np.nan)
            locus = ""
            if depth > 0:   # best 25-kb window of the profile (where a missing chromosomal piece belongs)
                k5 = 5
                n5 = len(r) // k5
                b = int(np.argmax(r[:n5 * k5].reshape(n5, k5).sum(1))) * k5
                locus = f"{P.names[P.chrom[b]]}:{int(P.bin_start[b]):,}-{int(P.bin_start[b]) + k5 * res:,}"
            cache[mol] = dict(r=r, depth=int(depth), p=p, intra=intra, eff=eff, ac=ac, feats=feats, locus=locus)
        m = cache[mol]
        flags = []
        if m["depth"] < c["min_depth"]: flags.append(f"depth<{c['min_depth']}")
        if m["intra"] < c["min_intra"]: flags.append("no intra contacts")
        if m["depth"] == 0:
            flags.append("no contact with the host genome (other organism / not nuclear?)")
        elif m["intra"] > c["max_intra"]:
            flags.append("organelle-like (intra>50%)")
        cis = m["feats"][NAMES.index("max_window_frac")]
        if m["depth"] >= c["min_depth"] and cis == cis and cis >= c["cis_frac"]:
            flags.append(f"cis-anchored: {cis:.0%} of contacts near {m['locus']}")
        if ctype.get(n) == "organelle": flags.append("annotated organelle")
        if m["eff"] == m["eff"] and m["eff"] < c["min_effect"]: flags.append(f"weak amplitude ({m['eff']:.2f})")
        if big and comm_frac >= c["community_frac"]:
            call = "CHROMOSOMAL_FRAGMENT" if comm_host[k] > 0 else "OTHER_GENOME"
        elif m["depth"] >= c["min_depth"] and cis == cis and cis >= c["cis_frac"]:
            call = "CHROMOSOMAL_FRAGMENT"
        elif m["p"] != m["p"]:
            call = "NOT_EVALUABLE"
        elif m["p"] <= thr:
            call = "OTHER"
        else:
            qc_ok = c["min_intra"] <= m["intra"] <= c["max_intra"] and ctype.get(n) != "organelle"
            call = ("STRONG_CANDIDATE" if m["eff"] >= c["min_effect"] else "WEAK_SIGNAL") if qc_ok else "QC_FAIL"
        meta = info.loc[n].to_dict() if n in info.index else {}
        rows.append(dict(
            id=n, call=call, p_hitchhiker=round(m["p"], 4) if m["p"] == m["p"] else np.nan,
            origin=meta.get("source") or corigin.get(n, ""), sequence_type=ctype.get(n, ""), length=int(P.lengths[i]),
            gc=meta.get("gc", ""), coverage=meta.get("coverage", ""), coverage_type=meta.get("coverage_type", ""),
            relative_copy_number=meta.get("relative_copy_number", ""), circular=meta.get("circular", ""),
            circular_evidence=meta.get("circular_evidence", ""), ref_cov_fraction=meta.get("ref_cov_fraction", ""),
            ref_best_hit=meta.get("ref_best_hit", ""), original_name=meta.get("original_name", ""),
            hic_trans_contacts=int(tot_i), molecule=",".join(mol), molecule_n_contigs=len(mol),
            molecule_trans_contacts=m["depth"], intra_fraction=round(m["intra"], 4),
            effect_size=round(m["eff"], 3) if m["eff"] == m["eff"] else np.nan,
            autocorr1=round(m["ac"], 3) if m["ac"] == m["ac"] else np.nan,
            community=int(k), community_bp=int(comm_len[k]), community_has_host=bool(comm_host[k] > 0),
            community_contact_share=round(float(comm_frac), 3), best_host_window=m["locus"], flags="; ".join(flags),
            **{f"f_{a}": round(float(b), 4) for a, b in zip(NAMES, m["feats"])}))

    order = {"STRONG_CANDIDATE": 0, "WEAK_SIGNAL": 1, "QC_FAIL": 2, "OTHER": 3, "CHROMOSOMAL_FRAGMENT": 4,
             "OTHER_GENOME": 5, "NOT_EVALUABLE": 6}
    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values(["call", "p_hitchhiker"], key=lambda s: s.map(order) if s.name == "call" else -s.fillna(-1))
    df.to_csv(lay.report, sep="\t", index=False)

    # one JSON file per contig, profile plots
    if c["contig_files"] and len(df):
        d = os.path.join(lay.class_dir, "contigs")
        os.makedirs(d, exist_ok=True)
        for rec in df.to_dict("records"):
            rec = {k: (None if isinstance(v, float) and v != v else v) for k, v in rec.items()}
            rec.update(unit=lay.unit, species=lay.row["name"], reference=lay.row["ref"], resolution=res,
                       model=model.spec.get("name"))
            with open(os.path.join(d, f"{rec['id']}.json"), "w") as fh:
                json.dump(rec, fh, indent=1, default=str)
    if c["plots"] != "none" and len(df):
        d = os.path.join(lay.class_dir, "plots")
        os.makedirs(d, exist_ok=True)
        sel = df if c["plots"] == "all" else df[df.call.isin(["STRONG_CANDIDATE", "WEAK_SIGNAL", "QC_FAIL"])]
        done = set()
        for rec in sel.to_dict("records"):
            mol = tuple(rec["molecule"].split(","))
            if mol in done or cache[mol]["depth"] == 0:
                continue
            done.add(mol)
            plot_profile(cache[mol]["r"], P, f"{lay.row['name']} | {rec['molecule']} | {rec['call']} | p={rec['p_hitchhiker']} "
                         f"| {cache[mol]['depth']} contacts", os.path.join(d, f"{rec['id']}.png"))
    summary = df.call.value_counts().to_dict() if len(df) else {}
    with open(os.path.join(lay.class_dir, "summary.json"), "w") as fh:
        json.dump(dict(unit=lay.unit, n_candidates=len(df), calls=summary,
                       strong=df[df.call == "STRONG_CANDIDATE"].molecule.unique().tolist() if len(df) else []), fh, indent=1)
    logger.info(f"[classify:{lay.unit}] {summary}")
    return summary
