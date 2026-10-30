"""SLURM submission: one job per (step, species) for 'load' and per (step, unit) for the other steps,
chained with --dependency=afterok so that the whole pipeline is submitted at once.

Every job simply runs the corresponding kontiguity command on its unit; parameters are read from
<outdir>/kontiguity_config.yaml, so a job reproduces exactly the configuration of the submission.
"""
import os
import shlex
import subprocess
import sys
import time

import pandas as pd

from kontiguity.utils.runner import StepError


def _fname(x):
    return x.replace("/", "__")


def write_script(outdir, step, target, cfg, extra_args=""):
    s = cfg["slurm"]
    res = s["resources"][step]
    d = os.path.join(outdir, "slurm")
    os.makedirs(os.path.join(d, "scripts"), exist_ok=True)
    os.makedirs(os.path.join(d, "logs"), exist_ok=True)
    name = f"{step}__{_fname(target)}"
    path = os.path.join(d, "scripts", name + ".sh")
    lines = ["#!/bin/bash", f"#SBATCH --job-name=kt_{step}_{_fname(target)[:40]}", f"#SBATCH --cpus-per-task={res['cpus']}",
             f"#SBATCH --mem={res['mem']}", f"#SBATCH --time={res['time']}",
             f"#SBATCH --output={d}/logs/{name}.%j.log", f"#SBATCH --error={d}/logs/{name}.%j.log"]
    if s.get("partition"): lines.append(f"#SBATCH --partition={s['partition']}")
    if s.get("qos"): lines.append(f"#SBATCH --qos={s['qos']}")
    if s.get("account"): lines.append(f"#SBATCH --account={s['account']}")
    for x in s.get("extra") or []:
        lines.append(f"#SBATCH {x}")
    lines += ["", "set -eo pipefail"] + list(s.get("setup") or [])
    py = shlex.quote(sys.executable)
    lines.append(f"{py} -m kontiguity {step} -o {shlex.quote(outdir)} --unit {shlex.quote(target)} "
                 f"-t ${{SLURM_CPUS_PER_TASK:-{res['cpus']}}} --no-slurm {extra_args}".rstrip())
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return path


def sbatch(script, deps=(), dry=False):
    cmd = ["sbatch", "--parsable"]
    deps = [d for d in deps if d]
    if deps:
        cmd.append("--dependency=afterok:" + ":".join(deps))
    cmd.append(script)
    if dry:
        return "DRY" + str(abs(hash(script)) % 100000)
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise StepError(f"sbatch failed for {script}: {p.stderr.strip()}")
    return p.stdout.strip().split(";")[0]


def record(outdir, rows, dry=False):
    path = os.path.join(outdir, "slurm", "jobs.dry_run.tsv" if dry else "jobs.tsv")
    df = pd.DataFrame(rows)
    if os.path.isfile(path):
        df = pd.concat([pd.read_csv(path, sep="\t", dtype=str), df], ignore_index=True)
    df.to_csv(path, sep="\t", index=False)
    return path


def wait(job_ids, poll=60, logger=None):
    ids = [j for j in job_ids if j and not j.startswith("DRY")]
    while ids:
        p = subprocess.run(["squeue", "-h", "-o", "%i", "-j", ",".join(ids)], capture_output=True, text=True)
        running = set(p.stdout.split()) if p.returncode == 0 else set()
        ids = [j for j in ids if j in running]
        if ids:
            if logger:
                logger.info(f"[slurm] {len(ids)} job(s) still queued/running")
            time.sleep(poll)
    states = {}
    p = subprocess.run(["sacct", "-n", "-X", "-P", "-o", "JobID,State", "-j", ",".join(job_ids)], capture_output=True, text=True)
    for line in p.stdout.splitlines():
        if "|" in line:
            j, st = line.split("|", 1)
            states[j] = st
    return states
