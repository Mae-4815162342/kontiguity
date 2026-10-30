"""Darwin Tree of Life data portal -> kontiguity sample table.

Keeps the species that have at least one assembly, one Hi-C run and one WGS run, and writes one row per
(species, assembly) with the DToL headers:
    name, ref, assembly_name, wgs, layout_WGS, protocol_WGS, platform_WGS, project_id_WGS,
    hic, layout_Hi-C, protocol_Hi-C, project_id_Hi-C, tax_id
Several runs are ';'-separated in one cell (the same table can be passed to every command with --table).
"""
import json
import re
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import requests

from kontiguity.samples import infer_tech

COLUMNS = ["name", "ref", "assembly_name", "wgs", "layout_WGS", "protocol_WGS", "platform_WGS", "project_id_WGS",
           "hic", "layout_Hi-C", "protocol_Hi-C", "project_id_Hi-C", "tax_id"]


def _get(url, params, timeout=60, tries=4):
    last = None
    for _ in range(tries):
        try:
            r = requests.get(url, params=params, timeout=timeout)
            r.raise_for_status()
            return r.json()
        except (requests.RequestException, ValueError) as e:
            last = e
    raise RuntimeError(f"DToL API request failed ({url}, {params}): {last}")


def fetch_all(url, page=100, threads=8, max_records=None, logger=None):
    first = _get(url, dict(limit=1, offset=0, sort="currentStatus:asc", current_class="kingdom"))
    count = first.get("count", 0)
    if max_records:
        count = min(count, max_records)
    offsets = list(range(0, count, page))
    if logger:
        logger.info(f"[dtol] {first.get('count', 0)} records on the portal, fetching {count}")
    with ThreadPoolExecutor(max_workers=threads) as ex:
        pages = list(ex.map(lambda o: _get(url, dict(limit=page, offset=o, sort="currentStatus:asc",
                                                     current_class="kingdom")).get("results", []), offsets))
    return [r for p in pages for r in p][:count]


def _txt(x):
    return "" if x is None else str(x).strip()


def records_to_table(records, ref_policy="all", primary_only=False, wgs_techs=None, species=None):
    rows = []
    pat = re.compile("|".join(species), re.I) if species else None
    for rec in records:
        src = rec.get("_source", rec)
        name = _txt(rec.get("_id") or src.get("organism") or src.get("scientific_name"))
        if pat and not pat.search(name):
            continue
        wgs, hic = [], []
        for e in src.get("experiment", []) or []:
            strat = _txt(e.get("library_strategy"))
            run = _txt(e.get("run_accession"))
            if not run:
                continue
            item = dict(run=run, layout=_txt(e.get("library_layout")),
                        protocol=_txt(e.get("library_construction_protocol")) or strat,
                        platform=" ".join(filter(None, [_txt(e.get("instrument_platform")), _txt(e.get("instrument_model"))])),
                        project=_txt(e.get("study_accession")))
            if strat == "Hi-C":
                hic.append(item)
            elif strat == "WGS":
                item["tech"] = infer_tech(item["protocol"], item["platform"], item["layout"])
                if not wgs_techs or item["tech"] in wgs_techs:
                    wgs.append(item)
        if not wgs or not hic:
            continue
        assemblies = []
        for a in src.get("assemblies", []) or []:
            acc = _txt(a.get("accession"))
            if acc:
                ver = _txt(a.get("version"))
                acc = acc if "." in acc or not ver else f"{acc}.{ver}"
                assemblies.append(dict(acc=acc, name=_txt(a.get("assembly_name") or a.get("description") or "")))
        annotated = {_txt(a.get("accession")) for a in (src.get("annotation", []) or [])}
        if ref_policy == "annotated":
            assemblies = [a for a in assemblies if a["acc"] in annotated or a["acc"].split(".")[0] in {x.split(".")[0] for x in annotated}]
            if not assemblies:
                assemblies = [dict(acc=a, name="") for a in sorted(annotated) if a]
        if primary_only:
            assemblies = [a for a in assemblies if not re.search(r"alternate|haplotype 2|hap2", a["name"], re.I)]
        if ref_policy == "latest" and assemblies:
            best = {}
            for a in assemblies:   # latest version of every assembly accession
                base, _, v = a["acc"].partition(".")
                if base not in best or (v or "0") > (best[base]["acc"].partition(".")[2] or "0"):
                    best[base] = a
            assemblies = [max(best.values(), key=lambda a: a["acc"])]
        seen = set()
        for a in assemblies:
            if a["acc"] in seen:
                continue
            seen.add(a["acc"])
            rows.append(dict(
                name=name, ref=a["acc"], assembly_name=a["name"],
                wgs=";".join(w["run"] for w in wgs), layout_WGS=";".join(w["layout"] for w in wgs),
                protocol_WGS=";".join(w["protocol"] for w in wgs), platform_WGS=";".join(w["platform"] for w in wgs),
                project_id_WGS=";".join(sorted({w["project"] for w in wgs})),
                hic=";".join(h["run"] for h in hic), **{"layout_Hi-C": ";".join(h["layout"] for h in hic),
                                                         "protocol_Hi-C": ";".join(h["protocol"] for h in hic),
                                                         "project_id_Hi-C": ";".join(sorted({h["project"] for h in hic}))},
                tax_id=_txt(src.get("tax_id") or src.get("taxId"))))
    return pd.DataFrame(rows, columns=COLUMNS)


def dtol_table(out, url, ref_policy="all", primary_only=False, wgs_techs=None, species=None, max_records=None,
               threads=8, json_dump=None, from_json=None, logger=None):
    if from_json:
        with open(from_json) as fh:
            data = json.load(fh)
        records = data.get("results", data) if isinstance(data, dict) else data
    else:
        records = fetch_all(url, threads=threads, max_records=max_records, logger=logger)
        if json_dump:
            with open(json_dump, "w") as fh:
                json.dump(records, fh)
    df = records_to_table(records, ref_policy, primary_only, wgs_techs, species)
    sep = "\t" if out.endswith((".tsv", ".txt")) else ","
    df.to_csv(out, sep=sep, index=False)
    if logger:
        logger.info(f"[dtol] {df.name.nunique()} species / {len(df)} (species, assembly) rows written to {out}")
    return df
