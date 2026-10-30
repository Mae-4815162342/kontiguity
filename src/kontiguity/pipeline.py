"""Execution engine shared by every command (local or SLURM)."""
import os
import traceback

import pandas as pd

from kontiguity import slurm
from kontiguity.settings import STEPS
from kontiguity.utils.runner import Layout, StepError, set_tmp_dir


def _step_fn(step):
    if step == "load":
        from kontiguity.load import load_species
        return load_species
    if step == "retrieve":
        from kontiguity.retrieve import retrieve_unit
        return retrieve_unit
    if step == "filter":
        from kontiguity.filter import filter_unit
        return filter_unit
    if step == "map":
        from kontiguity.hicmap import map_unit
        return map_unit
    if step == "classify":
        from kontiguity.classification.classify import classify_unit
        return classify_unit
    if step == "describe":
        from kontiguity.describe_step import describe_unit
        return describe_unit
    raise ValueError(step)


def targets(samples, step):
    """'load' works per species (reads are shared by all references); other steps per unit."""
    if step == "load":
        return list(dict.fromkeys(samples.species))
    return list(samples.unit)


# inputs a step needs from 'load' when it is run on its own
NEEDS = {"retrieve": ("ref", "wgs"), "filter": ("ref",), "map": ("hic",), "classify": (), "describe": ()}


def ensure_inputs(step, lay, samples, cfg, logger, user_inputs):
    """Prepares, if needed, what a step run on its own requires (reference, reads, user files)."""
    import shutil
    from kontiguity.load import load_species
    from kontiguity.load import missing_inputs
    need = set(NEEDS.get(step, ()))
    if step == "map" and not user_inputs.get("genome"):
        need |= {"ref"}
    if step == "filter" and user_inputs.get("contigs"):
        need = {"ref"}
    if need and not lay.is_done("load"):
        why = "'load' not done"
    elif need and missing_inputs(lay, need):   # e.g. reads deleted by --no-tmp, step re-run with --force
        need = missing_inputs(lay, need) | (need & {"ref"})
        why = "downloaded files missing (deleted by --no-tmp?)"
    else:
        need = set()
    if need:
        sp = [Layout(outdir_of(lay), r) for r in samples[samples.species == lay.species].to_dict("records")]
        logger.info(f"[{step}:{lay.unit}] {why}: preparing {', '.join(sorted(need))}")
        load_species(sp, cfg, logger, need=tuple(need))
    stage_user_files(step, lay, cfg, user_inputs)


def stage_user_files(step, lay, cfg, user_inputs):
    """Copies user-provided intermediate files where the step expects them."""
    import shutil
    if step == "filter" and user_inputs.get("contigs"):
        import pandas as pd
        from kontiguity.fasta import read_fasta, write_record, gc
        os.makedirs(lay.retrieve_dir, exist_ok=True)
        rows = []
        with open(lay.raw_contigs, "w") as out:
            for i, (h, s) in enumerate(read_fasta(user_inputs["contigs"])):
                cid = f"user_{i}"
                write_record(out, cid, s)
                rows.append(dict(id=cid, source="user", original_name=h.split()[0], length=len(s), gc=round(gc(s), 4),
                                 coverage=float("nan"), coverage_type="", circular="unknown", circular_evidence=""))
        pd.DataFrame(rows).to_csv(lay.raw_info, sep="\t", index=False)
    if step == "map" and user_inputs.get("genome"):
        os.makedirs(lay.contigs_dir, exist_ok=True)
        shutil.copy(user_inputs["genome"], lay.genome_fa)
        _types_table(lay.genome_fa, lay.genome_chroms, user_inputs.get("chroms"), cfg)
    if step == "classify" and user_inputs.get("mcool"):
        os.makedirs(lay.hic_dir, exist_ok=True)
        if os.path.lexists(lay.mcool):
            os.remove(lay.mcool)
        os.symlink(os.path.abspath(user_inputs["mcool"]), lay.mcool)
        if not os.path.isfile(lay.genome_chroms) or user_inputs.get("chroms"):
            _types_table(lay.mcool, lay.genome_chroms, user_inputs.get("chroms"), cfg)


def _types_table(source, out, chroms, cfg):
    """Sequence-type table of a user genome/mcool: from --chroms, otherwise sequences >= host_min_len are
    chromosomes and the others are candidate contigs."""
    import pandas as pd
    if source.endswith((".mcool", ".cool")):
        import h5py
        with h5py.File(source, "r") as f:
            g = f["resolutions"][sorted(f["resolutions"], key=int)[0]] if "resolutions" in f else f
            names = [x.decode() for x in g["chroms/name"][:]]
            lens = [int(x) for x in g["chroms/length"][:]]
    else:
        from kontiguity.fasta import lengths
        ln = lengths(source)
        names, lens = list(ln), list(ln.values())
    t = pd.DataFrame(dict(id=names, length=lens))
    user = {}
    if chroms:
        from kontiguity.load import read_chroms_table
        user = read_chroms_table(chroms)
    big = cfg["classify"]["host_min_len"]
    t["sequence_type"] = [user[n][0] if n in user else ("chromosome" if l >= big and not chroms else "recovered_contig")
                          for n, l in zip(t.id, t.length)]
    t["sequence_name"] = [user[n][1] if n in user else n for n in t.id]
    t["origin"] = ["reference" if n in user or l >= big else "user" for n, l in zip(t.id, t.length)]
    t["description"] = ""
    os.makedirs(os.path.dirname(out), exist_ok=True)
    t.to_csv(out, sep="\t", index=False)


def _rm(path):
    import shutil
    if os.path.islink(path) or os.path.isfile(path):
        os.remove(path)
    elif os.path.isdir(path):
        shutil.rmtree(path, ignore_errors=True)


def _du(paths):
    tot = 0
    for p in paths:
        if os.path.isfile(p) and not os.path.islink(p):
            tot += os.path.getsize(p)
        for root, _, files in os.walk(p):
            for f in files:
                fp = os.path.join(root, f)
                if not os.path.islink(fp):
                    try:
                        tot += os.path.getsize(fp)
                    except OSError:
                        pass
    return tot


def cleanup(step, lay, cfg, logger):
    """--no-tmp: deletes the bulky files of a unit as soon as no later step needs them.

    retrieve -> read bins, indexes, assemblies; WGS reads and Logan contigs of the species once every
                reference of the species is retrieved
    map      -> hicstuff alignments / pairs; Hi-C reads of the species once every reference is mapped
    describe -> the .mcool (last user of the contact map)
    A step re-run later (--force) downloads again what it needs."""
    if not cfg.get("no_tmp"):
        return
    import glob
    from kontiguity.samples import load_samples
    gone = []
    try:
        all_units = load_samples(lay.outdir)
        sibs = [Layout(lay.outdir, r) for r in all_units[all_units.species == lay.species].to_dict("records")]
    except Exception:   # no samples.tsv (unit tests): this unit only
        sibs = [lay]
    if step == "retrieve":
        gone += [lay.retrieve_work] if lay.tmp else (
            glob.glob(os.path.join(lay.retrieve_work, "unmapped")) + glob.glob(os.path.join(lay.retrieve_work, "*.mmi")) +
            glob.glob(os.path.join(lay.retrieve_work, "bt2")) + glob.glob(os.path.join(lay.retrieve_work, "assembly_*")))
        if all(s.is_done("retrieve") or s.unit == lay.unit for s in sibs):
            gone += [lay.wgs_dir, lay.logan_dir]
    elif step == "map":
        gone += [lay.hic_work]
        if all(s.is_done("map") or s.unit == lay.unit for s in sibs):
            gone += [lay.hic_dir_reads]
    elif step == "describe":
        if os.path.islink(lay.mcool) and os.path.realpath(lay.mcool) == os.path.realpath(lay.mcool_store):
            gone += [lay.mcool_store, lay.mcool]
        elif not os.path.islink(lay.mcool):   # a link to a user-provided map is never followed
            gone += [lay.mcool]
    gone = [g for g in gone if os.path.lexists(g)]
    if not gone:
        return
    size = _du(gone)
    for g in gone:
        _rm(g)
        if lay.tmp:   # empty parent folders of the temporary space
            d = os.path.dirname(g)
            while d.startswith(lay.tmp + os.sep):
                try:
                    os.rmdir(d)
                except OSError:
                    break
                d = os.path.dirname(d)
    human = f"{size / 1e9:.1f} GB" if size >= 1e8 else f"{size / 1e6:.1f} MB"
    logger.info(f"[{step}:{lay.unit}] --no-tmp: deleted {human} ({', '.join(gone)})")


def outdir_of(lay):
    return lay.outdir


def use_tmp_dir(outdir, cfg):
    """Registers --tmp-dir for this process; tools that honour TMPDIR (sort, cooler...) write there too."""
    root = set_tmp_dir(outdir, cfg.get("tmp_dir"))
    if root:
        os.makedirs(os.path.join(root, "TMPDIR"), exist_ok=True)
        os.environ["TMPDIR"] = os.path.join(root, "TMPDIR")
    return root


def run_local(step, samples, outdir, cfg, logger, force=False, user_inputs=None):
    use_tmp_dir(outdir, cfg)
    fn = _step_fn(step)
    failed = []
    for tgt in targets(samples, step):
        rows = samples[samples.species == tgt] if step == "load" else samples[samples.unit == tgt]
        lays = [Layout(outdir, r) for r in rows.to_dict("records")]
        lay = lays[0]
        if lay.is_done(step) and not force:
            logger.info(f"[{step}:{tgt}] already done (use --force to redo)")
            continue
        lay.mark(step, "running")
        try:
            if step != "load":
                ensure_inputs(step, lay, samples, cfg, logger, user_inputs or {})
            out = fn(lays, cfg, logger) if step == "load" else fn(lay, cfg, logger)
            lay.mark(step, "done", dict(result=str(out)))
            if step != "load":
                cleanup(step, lay, cfg, logger)
            logger.info(f"[{step}:{tgt}] done")
        except Exception as e:   # one failing unit must not stop the others
            msg = str(e) if isinstance(e, StepError) else traceback.format_exc()
            lay.mark(step, "failed", dict(error=msg[-4000:]))
            logger.error(f"[{step}:{tgt}] FAILED: {msg}")
            failed.append(tgt)
    return failed


def run_pipeline_local(samples, outdir, cfg, logger, first="load", last="describe", force=False, user_inputs=None):
    steps = STEPS[STEPS.index(first): STEPS.index(last) + 1]
    failed_units = set()
    for step in steps:
        if step == "load":
            sel = samples
        else:
            sel = samples[~samples.unit.isin(failed_units) & ~samples.species.isin(failed_units)]
        # user files only replace the input of the FIRST step that is run
        ui = (user_inputs or {}) if step == first else {}
        failed_units |= set(run_local(step, sel, outdir, cfg, logger, force, ui))
    return sorted(failed_units)


def submit(samples, outdir, cfg, logger, first="load", last="describe", force=False, dry=False, user_inputs=None):
    """Submits the requested steps of every unit, chained with afterok dependencies."""
    steps = STEPS[STEPS.index(first): STEPS.index(last) + 1]
    use_tmp_dir(outdir, cfg)
    if user_inputs and first != "load":
        for r in samples.to_dict("records"):
            stage_user_files(first, Layout(outdir, r), cfg, user_inputs)
    extra = "--force" if force else ""
    jobs, rows = {}, []
    for step in steps:
        for tgt in targets(samples, step):
            if step == "load":
                deps = []
            else:
                sp = samples[samples.unit == tgt].species.iloc[0]
                prev = STEPS[STEPS.index(step) - 1]
                deps = [jobs.get((prev, sp if prev == "load" else tgt))]
            lay = Layout(outdir, (samples[samples.species == tgt] if step == "load" else samples[samples.unit == tgt]).iloc[0])
            if lay.is_done(step) and not force:
                logger.info(f"[slurm] {step}:{tgt} already done, not submitted")
                continue
            script = slurm.write_script(outdir, step, tgt, cfg, extra)
            jid = slurm.sbatch(script, deps, dry=dry)
            jobs[(step, tgt)] = jid
            rows.append(dict(job_id=jid, step=step, target=tgt, depends_on=",".join(d for d in deps if d), script=script))
            logger.info(f"[slurm] {step}:{tgt} -> job {jid}" + (f" (after {','.join(d for d in deps if d)})" if any(deps) else ""))
    if rows:
        slurm.record(outdir, rows, dry)
    if cfg["slurm"]["wait"] and not dry:
        states = slurm.wait([r["job_id"] for r in rows], cfg["slurm"]["poll_seconds"], logger)
        bad = {j: s for j, s in states.items() if not s.startswith("COMPLETED")}
        if bad:
            logger.error(f"[slurm] jobs not completed: {bad}")
            return sorted(bad)
    return []


def status_table(samples, outdir):
    rows = []
    for r in samples.to_dict("records"):
        lay = Layout(outdir, r)
        d = dict(unit=r["unit"])
        for step in STEPS:
            st = "pending"
            for s in ("done", "failed", "running"):
                if os.path.isfile(lay.marker(step, s)):
                    st = s
                    break
            d[step] = st
        if os.path.isfile(lay.report):
            rep = pd.read_csv(lay.report, sep="\t")
            d["strong_candidates"] = int((rep.call == "STRONG_CANDIDATE").sum())
            d["weak_signal"] = int((rep.call == "WEAK_SIGNAL").sum())
        rows.append(d)
    return pd.DataFrame(rows)
