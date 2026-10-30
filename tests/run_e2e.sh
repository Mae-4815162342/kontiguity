#!/bin/bash
# End-to-end tests on simulated data (needs the conda environment of environment.yml).
#   bash tests/run_e2e.sh /path/to/workdir
set -euo pipefail
W=${1:-test_results}; mkdir -p "$W"
HERE=$(cd "$(dirname "$0")" && pwd)
python "$HERE/sim/make_sim.py" "$HERE/../test_data/load/S_cerevisiae.fna" "$W/sim"
S=$W/sim
# 1. local files, short reads, de novo assembly of the read bin
kontiguity run -o "$W/short" -n "Saccharomyces simulatus" -r "$S/reference.fa" --chroms "$S/chromosomes.csv" \
    --wgs "$S/wgs_R1.fq.gz:$S/wgs_R2.fq.gz" --hic "$S/hic_R1.fq.gz:$S/hic_R2.fq.gz" --hic-protocol "Arima v2" --no-logan -t 4 \
    --known "$S/truth.fa"
# where did the simulated plasmid / bacteria / missing chromosome piece end up?
cat "$W"/short/*/*/describe/known_hits.tsv
# 2. HiFi reads (Flye + tandem circles), up to the new reference
kontiguity run -o "$W/hifi" -n "Saccharomyces simulatus" -r "$S/reference.fa" --chroms "$S/chromosomes.csv" \
    --wgs "$S/hifi.fq.gz" --wgs-tech hifi --hic "$S/hic_R1.fq.gz:$S/hic_R2.fq.gz" --no-logan -t 4 --to filter
# 3. network code paths (ENA reference, Logan contigs, SRA mirror + fasterq-dump, ENA fallback) on a local mock server
python "$HERE/mock_network_run.py" "$S" "$W/mock"
# 4. SLURM orchestration with fake sbatch/squeue/sacct (jobs run synchronously)
mkdir -p "$W/fake_slurm_state"
PATH="$HERE/fake_slurm:$PATH" FAKE_SLURM_DIR="$W/fake_slurm_state" kontiguity run -o "$W/slurm" -n "Saccharomyces simulatus" \
    -r "$S/reference.fa" --chroms "$S/chromosomes.csv" --wgs "$S/wgs_R1.fq.gz:$S/wgs_R2.fq.gz" \
    --hic "$S/hic_R1.fq.gz:$S/hic_R2.fq.gz" --no-logan --slurm --slurm-wait
for d in short hifi mock slurm; do kontiguity status -o "$W/$d"; done
