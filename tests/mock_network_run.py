"""End-to-end test of the network code paths against a local mock server (no internet needed).

  * reference GCA_TEST000001.1  -> mock ENA fasta API (ENA-style headers)
  * WGS run SRR0000001          -> mock Logan bucket (contigs exist => the reads must NOT be downloaded)
  * Hi-C run ERR9000001         -> mock SRA open-data mirror (wget) + fake fasterq-dump on the local file
  * Hi-C run ERR9000002         -> absent from the SRA mirror -> ENA fastq fallback

    python tests/mock_network_run.py SIM_DIR OUTDIR
"""
import functools
import gzip
import http.server
import os
import shutil
import stat
import subprocess
import sys
import threading

SIM, OUT = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
SRV = OUT + "_server"
if not os.environ.get("MOCK_KEEP_OUT"):   # MOCK_KEEP_OUT=1: re-run on an existing output
    shutil.rmtree(OUT, ignore_errors=True)
shutil.rmtree(SRV, ignore_errors=True)
os.makedirs(SRV)

# --------------------------------------------------------------- served tree
os.makedirs(f"{SRV}/logan/c/SRR0000001")
shutil.copy(f"{SIM}/logan/c/SRR0000001/SRR0000001.contigs.fa.zst", f"{SRV}/logan/c/SRR0000001/")
os.makedirs(f"{SRV}/ena/fasta")
names = {"NC_001133.9": "chromosome: I", "NC_001134.8": "chromosome: II", "NC_001135.5": "chromosome: III",
         "NC_001138.5": "chromosome: VI", "NC_001224.1": "mitochondrion"}
with open(f"{SIM}/reference.fa") as fh, gzip.open(f"{SRV}/ena/fasta/GCA_TEST000001.1.fa.gz", "wt") as out:
    for line in fh:
        if line.startswith(">"):
            acc = line[1:].split()[0]
            out.write(f">ENA|{acc.split('.')[0]}|{acc} Saccharomyces simulatus genome assembly, {names[acc]}\n")
        else:
            out.write(line)
# SRA object for ERR9000001 (content is irrelevant: the fake fasterq-dump reads fixtures)
os.makedirs(f"{SRV}/sra/ERR9000001")
open(f"{SRV}/sra/ERR9000001/ERR9000001", "w").write("fake sra object\n")
# ENA fastq for ERR9000002 = second half of the Hi-C reads
os.makedirs(f"{SRV}/ena/files")
half = {}
for r in (1, 2):
    with gzip.open(f"{SIM}/hic_R{r}.fq.gz", "rt") as fh:
        lines = fh.readlines()
    n = len(lines) // 8 * 4
    half[r] = (lines[:n], lines[n:])
    with gzip.open(f"{SRV}/ena/files/ERR9000002_{r}.fastq.gz", "wt") as out:
        out.writelines(half[r][1])

# --------------------------------------------------------------- server
handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=SRV)
httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
port = httpd.server_address[1]
threading.Thread(target=httpd.serve_forever, daemon=True).start()
base = f"http://127.0.0.1:{port}"
with open(f"{SRV}/ena/filereport_ERR9000002.tsv", "w") as fh:
    fh.write(f"run_accession\tfastq_ftp\nERR9000002\t{base}/ena/files/ERR9000002_1.fastq.gz;{base}/ena/files/ERR9000002_2.fastq.gz\n")
with open(f"{SRV}/ena/filereport_ERR9000001.tsv", "w") as fh:
    fh.write("run_accession\tfastq_ftp\nERR9000001\t\n")

# --------------------------------------------------------------- fake fasterq-dump (works on the local file only)
fixtures = f"{SRV}/fixtures"
os.makedirs(fixtures)
for r in (1, 2):
    with open(f"{fixtures}/R{r}.fastq", "w") as out:
        out.writelines(half[r][0])
bindir = f"{SRV}/bin"
os.makedirs(bindir)
fake = f"""#!/bin/bash
# fake fasterq-dump: checks it was given a LOCAL file (no prefetch / no network), writes <acc>_1/_2.fastq
out=.; last=""
while [ $# -gt 0 ]; do case "$1" in -O) out=$2; shift 2;; -e|--temp) shift 2;; -*) shift;; *) last=$1; shift;; esac; done
[ -f "$last" ] || {{ echo "not a local file: $last" >&2; exit 3; }}
acc=$(basename "$last")
cp {fixtures}/R1.fastq "$out/${{acc}}_1.fastq"; cp {fixtures}/R2.fastq "$out/${{acc}}_2.fastq"
echo "spots read: fake" >&2
"""
open(f"{bindir}/fasterq-dump", "w").write(fake)
os.chmod(f"{bindir}/fasterq-dump", os.stat(f"{bindir}/fasterq-dump").st_mode | stat.S_IEXEC)

cfg = f"""urls:
  sra: "{base}/sra/{{acc}}/{{acc}}"
  ena_filereport: "{base}/ena/filereport_{{acc}}.tsv"
  ena_fasta: "{base}/ena/fasta/{{acc}}.fa.gz"
logan:
  url: "{base}/logan/c/{{acc}}/{{acc}}.contigs.fa.zst"
"""
os.makedirs(OUT, exist_ok=True)
open(f"{OUT}.yaml", "w").write(cfg)
table = f"{OUT}.table.tsv"
open(table, "w").write("name\tref\twgs\tlayout_WGS\tprotocol_WGS\thic\tlayout_Hi-C\tprotocol_Hi-C\n"
                       "Saccharomyces simulatus\tGCA_TEST000001.1\tSRR0000001\tPAIRED\tIllumina NovaSeq\t"
                       "ERR9000001;ERR9000002\tPAIRED;PAIRED\tHi-C - Arima v2\n")
env = dict(os.environ, PATH=bindir + os.pathsep + os.environ["PATH"], NO_PROXY="127.0.0.1,localhost",
           no_proxy="127.0.0.1,localhost")
cmd = ["kontiguity", "run", "-o", OUT, "--table", table, "--config", f"{OUT}.yaml", "-t", "2"] + sys.argv[3:]
print(" ".join(cmd), flush=True)
rc = subprocess.call(cmd, env=env)
httpd.shutdown()
sys.exit(rc)
