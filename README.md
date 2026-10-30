# Kontiguity

Kontiguity finds **extrachromosomal DNA molecules in eukaryotic genomic datasets** (plasmids first of all) and
classifies them from their **Hi-C contact signal**.

For every (species, reference assembly) pair it:

1. gathers candidate contigs, either from **Logan** (assemblies of the whole SRA) or, when Logan does not
   have the run, by assembling the WGS reads that do **not** align on the reference (the "read bin");
2. re-aligns these contigs on the reference and keeps those that do not match it;
3. builds a **new reference** (reference + kept contigs);
4. maps the **Hi-C** reads on this new reference (hicstuff + cooler);
5. **classifies** the contact profile of every contig with a model trained on the natural eukaryotic
   plasmids whose behaviour is known (2-micron of *S. cerevisiae*, its *Lachancea* relative, Ddp5 of
   *D. discoideum*: they hitchhike on silent chromosome regions, Girard et al. 2025);
6. **describes** the recovered contigs: one table of all of them, the search of known sequences (`--known`,
   e.g. "is the 2-micron plasmid among my contigs?") and the kontiguity descriptive Hi-C figures.

Datasets can come from the **Darwin Tree of Life** portal (species with an assembly, Hi-C and WGS runs), from
a table using the DToL headers, or from accessions / local files given on the command line. Every step is a
separate command, `kontiguity run` chains them, and everything can be submitted to **SLURM**.

---

## Contents
1. [Installation](#installation)
2. [Quick start](#quick-start)
3. [Input formats](#input-formats) ← read this first
4. [Commands](#commands)
5. [Output tree](#output-tree)
6. [The per-contig report](#the-per-contig-report)
7. [Logan](#logan)
8. [Downloads and SLURM clusters](#downloads-and-slurm-clusters)
9. [Configuration and defaults](#configuration-and-defaults)
10. [The classifier](#the-classifier)
11. [Tests](#tests)
12. [Changes from 0.0.x](#changes-from-00x)

---

## Installation

Kontiguity drives external programs. The simplest is the provided conda environment:

```bash
git clone https://github.com/Mae-4815162342/kontiguity.git
cd kontiguity
conda env create -f environment.yml        # minimap2, samtools, bowtie2, megahit, flye, sra-tools,
conda activate kontiguity                  # hicstuff, cooler, wget, zstd, pigz + python packages
pip install -e .                           # or: pip install kontiguity
kontiguity --help
```

| Program | Used by | For |
|---|---|---|
| `wget`, `fasterq-dump` (sra-tools ≥ 3), `zstd`, `pigz` (optional) | load | reads, Logan contigs |
| `minimap2`, `samtools` | retrieve, filter | read bin, contig re-alignment |
| `bowtie2` | retrieve (option), map | short-read aligner, hicstuff |
| `megahit` | retrieve | short-read assembly |
| `flye` ≥ 2.9 | retrieve | long-read assembly (HiFi, CLR, ONT) |
| `hicstuff` ≥ 3.2.5, `cooler` | map | Hi-C pairs, .mcool |

A step only checks the programs it needs, and fails with a clear message if one is missing.

---

## Quick start

```bash
# 1. Darwin Tree of Life: build the sample table, then run everything on SLURM
kontiguity dtol -o dtol.tsv --primary-only
kontiguity run -o results --table dtol.tsv --slurm --slurm-partition common \
    --slurm-setup "source activate kontiguity"

# 2. your own table (same headers, see below), locally with 16 threads
kontiguity run -o results --table my_samples.tsv -t 16

# 3. accessions directly
kontiguity run -o results -n "Vespula vulgaris" -r GCA_905475345.1 \
    --wgs ERR6635597 --hic ERR6054674,ERR6054675 --hic-protocol "Arima v2"

# 4. local files
kontiguity run -o results -n "My species" -r genome.fa --chroms genome_chromosomes.csv \
    --wgs reads_R1.fq.gz:reads_R2.fq.gz --hic hic_R1.fq.gz:hic_R2.fq.gz

# follow the progress, then read the reports
kontiguity status -o results
cat results/*/*/classification/contigs_report.tsv
```

---

## Input formats

Every command (`load`, `retrieve`, `filter`, `map`, `classify`, `run`) accepts the **same input options**.
They are given once: the normalized sample table is saved as `<outdir>/samples.tsv`, and later commands on
the same `--outdir` reuse it without any input option. Use `--unit` to restrict a command to some units.

**Unit.** The analysis unit is one *(species, reference assembly)* pair, named `Species_name/REF`
(e.g. `Vespula_vulgaris/GCA_905475345.1`). Reads and Logan contigs belong to the species and are shared by
all its references; one species with three assemblies gives three units and three new references.

### 1. Direct identifiers

| Option | Value | Example |
|---|---|---|
| `-n, --name` | species / experiment name (spaces allowed, replaced by `_` in paths) | `-n "Vespula vulgaris"` |
| `-r, --ref` | assembly accession(s) fetched from ENA **or** FASTA path(s); several comma-separated (one unit each) | `-r GCA_905475345.1,GCA_905475355.1` |
| `--wgs` | WGS ("assembly") runs: SRA/ENA run accessions (`SRR…`, `ERR…`, `DRR…`) **or** local FASTQ. Runs are comma-separated; the two files of a paired local run are joined by `:` | `--wgs ERR6688535,ERR6939243`<br>`--wgs a_R1.fq.gz:a_R2.fq.gz,long.fq.gz` |
| `--hic` | Hi-C runs, same syntax (Hi-C must be paired) | `--hic ERR6054674` or `--hic h_R1.fq.gz:h_R2.fq.gz` |
| `--wgs-tech` | technology of the `--wgs` runs: `hifi`, `clr`, `ont`, `short` or `auto`; one value per run, or one for all | `--wgs-tech short,hifi` |
| `--hic-protocol` | Hi-C library protocol, used to choose the restriction enzymes | `--hic-protocol "Arima v2"` |
| `--chroms` | sequence table of a **local** reference (see 3) | `--chroms chromosomes.csv` |

`auto` technology: paired runs are short reads; for local single-end files the mean read length decides
(> 1 kb: HiFi). For DToL/ENA runs the technology comes from the protocol/platform columns.

### 2. Sample table (`--table`, Darwin Tree of Life headers)

CSV or TSV (`.tsv`/`.txt` → tab-separated). This is exactly what `kontiguity dtol` writes.

| Column | Required | Content |
|---|---|---|
| `name` | **yes** | species name |
| `ref` | **yes** | assembly accession or FASTA path (several: `;`) |
| `wgs` | **yes** for retrieve | WGS run accession(s) or local FASTQ, several separated by `;` (paired local files: `R1:R2`) |
| `hic` | **yes** for map | Hi-C run accession(s) or local paired FASTQ (`R1:R2`), several separated by `;` |
| `layout_WGS` | no | `PAIRED` / `SINGLE`, one per WGS run (`;`) |
| `protocol_WGS` | no | e.g. `PacBio - HiFi`, `Chromium genome`, one per run: gives the technology |
| `platform_WGS` | no | e.g. `PACBIO_SMRT Sequel II`, `ILLUMINA`: also used for the technology |
| `wgs_tech` | no | explicit technology per run (`hifi;short`), overrides protocol/platform |
| `layout_Hi-C`, `protocol_Hi-C` | no | e.g. `PAIRED`, `Hi-C - Arima v2`: gives the enzymes |
| `enzymes` | no | explicit hicstuff enzymes for the row (`DpnII,HinfI`, `dnase`, …) |
| `chroms` | no | sequence table of a local reference (see 3) |
| `project_id_WGS`, `project_id_Hi-C`, `assembly_name`, `tax_id` | no | informative only |

Two layouts are accepted:

* **one row per (species, assembly)** with `;`-separated runs (what `kontiguity dtol` writes):
  ```
  name	ref	wgs	layout_WGS	protocol_WGS	hic	layout_Hi-C	protocol_Hi-C
  Hydraecia micacea	GCA_914767645.1	ERR6688535;ERR6939243	PAIRED;SINGLE	Chromium genome;PacBio - HiFi	ERR6688534	PAIRED	Hi-C - Arima v2
  ```
* **one row per WGS × Hi-C combination** (the format of kontiguity 0.0.x, e.g. `test_data/load/test_dataset.csv`):
  rows sharing `name` and `ref` are merged (union of the runs). Older tables where `wgs`/`hic` hold study
  accessions and the runs are in `run_acession_WGS`/`run_acession_Hi-C` are also understood.

Unknown columns are ignored, so a table exported from the DToL portal can be used after renaming its
columns.

### 3. Sequence table of a local reference (`--chroms` or column `chroms`)

Assemblies downloaded from ENA are typed from their FASTA headers
(`>ENA|OX637620|OX637620.1 … genome assembly, chromosome: 1`). For a **local** FASTA, give a CSV/TSV:

```
id,sequence_type,sequence_name
NC_001133.9,chromosome,I
NC_001134.8,chromosome,II
NC_001224.1,organelle,mitochondrion
```

`id` is the first word of the FASTA header. `sequence_type` follows ENA: `chromosome`, `organelle`,
`unlocalized`, `unplaced`, `scaffold`, `contig`, `plasmid`. **Chromosomes are the host** against which
contact profiles are computed; every other sequence is a candidate. Without a table, sequences
≥ 100 kb (`--host-min-len`) are taken as chromosomes.

### 4. Your own intermediate files

| Option | Command | Replaces |
|---|---|---|
| `--contigs contigs.fa` | `filter` (or `run --from filter`) | the output of `retrieve`: your own assembly / contigs are sorted against the reference |
| `--genome genome.fa` (+ `--chroms`) | `map` (or `run --from map`) | the new reference built by `filter` |
| `--mcool map.mcool` (+ `--chroms`) | `classify` | the Hi-C map built by `map`; contigs = non-chromosome sequences |

```bash
kontiguity classify -o out -n "S. cerevisiae" --mcool S_cerevisiae_contigs.mcool
```

A step launched on its own prepares the inputs it needs if `load` has not run (reference formatting,
reads download), so each command can be used independently.

---

## Commands

| Command | Scope (one job per) | Does |
|---|---|---|
| `kontiguity dtol` | – | queries the DToL portal and writes the sample table |
| `kontiguity load` | species | reference(s), Logan contigs, reads (only those needed) |
| `kontiguity retrieve` | unit | Logan contigs, or read bin → de novo contigs |
| `kontiguity filter` | unit | re-alignment on the reference, sorting, redundancy, circularity, new reference |
| `kontiguity map` | unit | hicstuff + cooler → `.mcool` on the new reference |
| `kontiguity classify` | unit | contact-profile classification → per-contig report |
| `kontiguity describe` | unit | table of the recovered contigs, known sequences (`--known`), descriptive figures |
| `kontiguity run` | all | `load → retrieve → filter → map → classify → describe` (`--from`, `--to`) |
| `kontiguity status` | – | done / running / failed / pending per unit and step, number of candidates |

Common options of every step and of `run`:

```
-o, --outdir DIR        output directory (required)
input options           --table | -n/-r/--wgs/--hic/--wgs-tech/--hic-protocol/--chroms   (see above)
-u, --unit TEXT         restrict to unit(s) / species (repeatable)
-t, --threads INT       threads per job [8]
--config FILE           YAML overriding any parameter (see Configuration)
--keep-tmp              keep intermediate files (read bin, assemblies, hicstuff tmp)
--tmp-dir DIR           put the bulky files (reads, Logan contigs, read bins, assemblies, Hi-C alignments,
                        .mcool) in DIR instead of the output tree; '$SCRATCH/x' allowed (see "Disk space")
--no-tmp                delete those files as soon as no later step needs them (.mcool after describe)
--force                 redo steps already done (steps are skipped when done, so runs can be resumed)
-v, --verbose           print the log on the terminal too (it is always in <outdir>/logs/kontiguity.log)
--slurm / SLURM options  see "Downloads and SLURM clusters"
```

Step options (all also accepted by `run`):

```
load      --no-logan  --logan-keep-reads  --logan-url  --download-source [auto|sra|ena]
          --keep-sra  --parallel-downloads
retrieve  --aligner [minimap2|bowtie2]  --min-contig-len  --megahit-opts "..."  --flye-opts "..."
filter    --contigs FILE  --max-ref-cov  --min-identity  --max-contig-len  --no-dedup  --minimap2-preset
map       --genome FILE  --enzymes  --hic-aligner [bowtie2|minimap2|bwa]  --hic-mapping [normal|iterative|cutsite]
          --hic-filter  --no-duplicates  --min-mapq  --resolutions 5000,10000,...  --no-balance
classify  --mcool FILE  --model  --threshold  --min-depth  --min-effect  --resolution  --host-min-len
          --recovered-only  --no-contig-files  --plots [none|candidates|all]
describe  --known FASTA  --describe-contigs [recovered|candidates|all]  --describe-max N
          --describe-binnings 10000,...  --formats png,pdf  --no-describe-plots
dtol      -o TABLE  --species REGEX  --ref-policy [all|latest|annotated]  --primary-only
          --wgs-tech hifi,short,...  --max-records N  --save-json F  --from-json F
```

`kontiguity dtol` keeps the species having at least one assembly, one Hi-C run and one WGS run, and writes
one row per (species, assembly): by default **every assembly** of a species becomes a reference
(`--ref-policy latest|annotated` to restrict, `--primary-only` to drop alternate haplotypes). `--save-json`
keeps the raw portal answer and `--from-json` rebuilds a table from it offline (e.g. on a cluster without
internet access).

---

## Output tree

```
<outdir>/
├── samples.tsv                     normalized sample table (one row per unit)
├── kontiguity_config.yaml          effective parameters (reused by later commands and SLURM jobs)
├── logs/kontiguity.log             log of every command
├── slurm/                          scripts/, logs/, jobs.tsv (job ids and dependencies)
└── <Species_name>/
    ├── reads/
    │   ├── manifest.tsv            every run: files, layout, technology, source (sra|ena|local|logan_only)
    │   ├── wgs/<RUN>_1.fastq.gz …  downloaded WGS reads (only runs without Logan contigs)
    │   └── hic/<RUN>_1.fastq.gz …  downloaded Hi-C reads
    ├── logan/<RUN>.contigs.fa      Logan contigs of each WGS run
    │                               (reads/wgs, reads/hic, logan/, retrieve working files, hic/hicstuff/ and
    │                                the .mcool go to --tmp-dir when given, see "Disk space")
    ├── logs/load.log, .status/load.done
    └── <REF>/                      one directory per reference = one unit
        ├── reference/
        │   ├── ref.fa              reference (sequences renamed to their accession)
        │   └── chromosomes.tsv     id, sequence_type, sequence_name, length, origin, description
        ├── retrieve/
        │   ├── contigs_raw.fa      all candidate contigs (Logan + de novo + tandem circles)
        │   ├── contigs_raw.tsv     their source, coverage and circularity
        │   ├── stats.json          host sequencing depth, runs used, assemblies
        │   └── assembly_<tech>/    MEGAHIT / Flye outputs (final contigs + logs)
        ├── contigs/
        │   ├── contigs.fa          kept contigs, renamed ktg_00001, ktg_00002, … (longest first)
        │   ├── contigs_info.tsv    one row per kept contig (see below)
        │   ├── contigs_vs_reference.paf
        │   ├── genome.fa           ← the NEW REFERENCE: reference + kept contigs
        │   └── genome_chromosomes.tsv  types of genome.fa (contigs: 'recovered_contig')
        ├── hic/
        │   ├── <Species>.<REF>.mcool   balanced, resolutions 5, 10, 50, 100 kb
        │   ├── chrom.sizes, hicstuff log, distance law
        │   └── hicstuff/           pairs (valid_idx_pcrfree.pairs.gz)
        ├── classification/
        │   ├── contigs_report.tsv  ← ONE ROW PER CONTIG, all the information below
        │   ├── contigs/<id>.json   one file per contig (same fields + unit, model, resolution)
        │   ├── plots/<id>.png      contact profile vs null envelope (candidates; --plots all for every contig)
        │   └── summary.json        number of contigs per call, strong candidates
        ├── describe/
        │   ├── recovered_contigs.tsv  every recovered contig: origin, length, GC, circularity, coverage, Hi-C, call
        │   ├── known_hits.tsv      (--known) matches of the known sequences in the reference and the contigs
        │   └── plots/              descriptive figures: <contig>.<bin>.summary.png, Mini_matrices.png, signals.csv
        ├── logs/<step>.log         commands and tool outputs of every step
        └── .status/<step>.done|failed|running
```

---

## The per-contig report

`classification/contigs_report.tsv` (and `classification/contigs/<id>.json`) gives, for every candidate:

| Field | Meaning |
|---|---|
| `call` | **STRONG_CANDIDATE**, WEAK_SIGNAL, QC_FAIL, OTHER, CHROMOSOMAL_FRAGMENT, OTHER_GENOME, NOT_EVALUABLE (see below) |
| `p_hitchhiker` | classifier probability of a 2-micron / Ddp5-like contact profile (molecule level) |
| `origin` | `logan_<RUN>`, `megahit_short`, `flye_hifi`, `readcircle_hifi`, `user`, or `reference` (scaffolds/organelles of the assembly) |
| `sequence_type`, `length`, `gc` | |
| `coverage`, `coverage_type` | Logan mean k-mer abundance (`ka`), Flye read depth, MEGAHIT k-mer multiplicity, or copies seen in reads |
| `relative_copy_number` | contig coverage / host coverage of the same source (Logan: median `ka` of the contigs matching the genome; reads: mapped bases / genome size). Approximate. |
| `circular`, `circular_evidence` | Logan graph self-loop, Flye circular contig, tandem copies in HiFi reads, or end self-overlap (trimmed) |
| `ref_cov_fraction`, `ref_best_hit` | part of the contig aligned on the reference, best reference sequence |
| `hic_trans_contacts` | Hi-C contacts of the contig with the host chromosomes |
| `molecule`, `molecule_n_contigs`, `molecule_trans_contacts` | contigs merged as one molecule (Hi-C community) and its contacts |
| `intra_fraction` | share of intra-molecule contacts (all real molecules: 2–35 %) |
| `effect_size` | amplitude of the signal beyond Poisson noise (known plasmids 0.35–0.71) |
| `autocorr1` | spatial autocorrelation of the profile (real hotspots span several bins) |
| `community`, `community_bp`, `community_has_host`, `community_contact_share` | Louvain community of the contig |
| `best_host_window` | host locus most contacted (for fragments: where they belong) |
| `flags` | QC messages |
| `f_*` | the 12 features given to the classifier |

**Calls**

| Call | Rule |
|---|---|
| `STRONG_CANDIDATE` | p > 0.5, amplitude ≥ 0.30, 0.5 % ≤ intra-contacts ≤ 50 %, not an annotated organelle |
| `WEAK_SIGNAL` | p > 0.5 but amplitude < 0.30 (e.g. organelles, rDNA episome of *E. histolytica*) |
| `QC_FAIL` | p > 0.5 but no intra-contig contacts, organelle-like, or annotated organelle |
| `OTHER` | p ≤ 0.5 |
| `CHROMOSOMAL_FRAGMENT` | ≥ 10 % of contacts in one 25-kb host window (missing piece of a chromosome), or contig inside a Mb-scale Hi-C community holding host chromosomes |
| `OTHER_GENOME` | contig inside a Mb-scale community without host chromosome (endosymbiont, contaminant) |
| `NOT_EVALUABLE` | < 500 trans contacts (flag *no contact with the host genome* when there is none at all) |

Fragments of one molecule (a plasmid split into several contigs) are grouped by Hi-C community detection and
scored together: e.g. the three Ddp5 contigs of *D. discoideum* form one molecule.

---

### Checking that a known molecule was recovered

```bash
kontiguity describe -o results --known 2micron.fa          # after a run, or as the last step of 'run'
cat results/Saccharomyces_cerevisiae/*/describe/known_hits.tsv
```
`known_hits.tsv` gives, for every known sequence, the reference sequences and recovered contigs it matches
(`known_covered`: share of the known sequence covered by that hit, `identity`, and the `call` /
`p_hitchhiker` of the contig). The log also says, for each known sequence, whether it was found among the
recovered contigs. Contigs matching a known sequence are drawn first in `describe/plots/`.

## Logan

[Logan](https://github.com/IndexThePlanet/Logan) provides assemblies (unitigs and contigs) of nearly all
public SRA runs. For each WGS run with an SRA/ENA accession, `load` checks
`https://s3.amazonaws.com/logan-pub/c/<RUN>/<RUN>.contigs.fa.zst` (HTTP, no AWS client needed):

* **contigs exist** → they are downloaded (wget) and decompressed; the reads of that run are **not
  downloaded** and no alignment/assembly is done: only the final re-alignment on each reference (`filter`)
  is performed. Their headers give the coverage (`ka:f:` mean k-mer abundance) and the assembly-graph links
  (`L:+:<n>:+`): a contig linked to itself closes a cycle and is reported **circular**.
* **contigs absent** → reads are downloaded and assembled.

`--no-logan` disables Logan (always assemble the reads); `--logan-keep-reads` uses both Logan contigs and the
reads; `--logan-url` changes the address (mirror).

---

## Downloads and SLURM clusters

**Reads.** `prefetch` does not work on most compute nodes, so it is never used. Reads are fetched with
`wget` from the NCBI SRA open-data mirror (`https://sra-pub-run-odp.s3.amazonaws.com/sra/<RUN>/<RUN>`), then
converted **locally** with `fasterq-dump` on the downloaded file, and gzipped. If the run is not on the
mirror, the ENA fastq.gz files are downloaded instead (`--download-source auto`, the default; `sra` or `ena`
to force one). Downloads resume (`wget -c`) and are retried. References come from the ENA browser API.

**SLURM.** Add `--slurm` to any command. `run --slurm` submits the whole pipeline at once:

* one `load` job per **species** (downloads are not duplicated across references),
* then, per unit, `retrieve → filter → map → classify → describe`, chained with `--dependency=afterok`.

Each job runs `python -m kontiguity <step> -o <outdir> --unit <unit>` and reads the parameters saved in
`<outdir>/kontiguity_config.yaml`, so jobs reproduce exactly the submitted configuration. Steps already done
are not resubmitted, failed units can be resubmitted with the same command.

```bash
kontiguity run -o results --table dtol.tsv --slurm \
    --slurm-partition common --slurm-qos normal --slurm-account myproject \
    --slurm-res retrieve=32,128G,72:00:00 --slurm-res map=32,96G,48:00:00 \
    --slurm-setup "module load conda" --slurm-setup "conda activate kontiguity" \
    --slurm-extra "--constraint=avx2"
kontiguity status -o results
```

`--dry-run` writes the scripts (`<outdir>/slurm/scripts/`) without submitting; `--slurm-wait` waits for the
jobs (polling `squeue`) and reports failures (`sacct`). Default resources per step:

| step | cpus | mem | time |
|---|---|---|---|
| load | 4 | 8G | 24 h |
| retrieve | 16 | 64G | 48 h |
| filter | 8 | 32G | 12 h |
| map | 16 | 48G | 48 h |
| classify | 4 | 16G | 6 h |
| describe | 2 | 16G | 6 h |

### Disk space: `--tmp-dir` and `--no-tmp`

A HiFi WGS run of a DToL species is typically 20–100 GB of reads and `fasterq-dump` needs about twice that
temporarily; Hi-C runs are of the same order; Logan contigs of a run are a few hundred MB to a few GB.
Runs covered by Logan avoid the WGS read download entirely. Without options, small intermediates (read-bin
indexes, assembly internals, hicstuff temporary files, `.sra`) are deleted unless `--keep-tmp` /
`--keep-sra`, but **downloaded reads, hicstuff alignments/pairs and the `.mcool` are kept**.

When the output space is small, send the bulky files elsewhere and delete them as you go:

```bash
kontiguity run --table dtol.tsv -o results --tmp-dir '$SCRATCH/kontiguity' --no-tmp --slurm ...
```

`--tmp-dir DIR` writes to `DIR/kontiguity_<outdir name>_<hash>/<species>/…`:

| moved to `--tmp-dir` | deleted by `--no-tmp` |
|---|---|
| `reads/wgs/`, `logan/` (downloaded WGS reads, Logan contigs, fasterq-dump temporary files) | once **every reference of the species** is retrieved |
| `retrieve/` working files (read bin, indexes, MEGAHIT/Flye assemblies) | after `retrieve` |
| `reads/hic/` (downloaded Hi-C reads) | once every reference of the species is mapped |
| `hic/hicstuff/` (BAM, pairs, pooled reads), base `.cool` | after `map` |
| the `.mcool` (`hic/<Species>.<REF>.mcool` in the output tree is then a link to it) | after `describe` |
| `$TMPDIR` of the jobs (cooler, sort…) | – |

* What stays in the output directory: references, `contigs_raw.fa`, `contigs/` (incl. the new reference
  `genome.fa`), hicstuff logs and distance law, `classification/`, `describe/`, logs. With `--no-tmp`, the
  output of a DToL species is a few MB plus its reference.
* Your **own** files (local fastq, `--mcool`, `--genome`, `--contigs`) are never moved nor deleted.
* Use quotes around `$VAR`: it is expanded inside each job, so the directory **must be visible from every
  node** (a shared scratch / work space). A per-job local disk (a `$TMPDIR` that is wiped when the job
  ends) does not work, because consecutive steps run in different jobs.
* Both options are saved in `kontiguity_config.yaml` like every parameter, so `status`, re-runs and SLURM
  jobs use the same locations. Do not change `--tmp-dir` in the middle of a run.
* Re-running a step after its inputs were deleted works: missing reads are downloaded again
  (`kontiguity run --from map --force …` re-downloads the Hi-C reads). `classify`/`describe` need the
  `.mcool`: with `--no-tmp`, re-run them with `--from map`.
* If the scratch space purges old files, a step whose files were purged is simply re-run the same way.

Compute nodes need outbound HTTPS for `load` (wget). If they have none, run `kontiguity load` on a login /
transfer node first, then submit the rest with `kontiguity run --from retrieve --slurm`.

---

## Configuration and defaults

All parameters can be set in a YAML file passed with `--config` (or with the command-line options). The
effective configuration is saved in `<outdir>/kontiguity_config.yaml`.

```yaml
threads: 16
logan:
  enabled: true
load:
  download_source: auto
retrieve:
  aligner: minimap2
  flye_opts: "--min-overlap 1000"
filter:
  max_ref_cov: 0.5
map:
  resolutions: [5000, 10000, 50000, 100000]
classify:
  min_depth: 500
slurm:
  partition: common
  setup: ["source activate kontiguity"]
  resources:
    map: {cpus: 32, mem: 96G, time: "72:00:00"}
urls:                       # mirrors
  sra: "https://sra-pub-run-odp.s3.amazonaws.com/sra/{acc}/{acc}"
```

Main defaults and why:

| Parameter | Default | Reason |
|---|---|---|
| read aligner | minimap2 (`sr`, `map-hifi`, `map-pb`, `map-ont` by technology) | one fast aligner for every read type |
| read bin (long reads) | unmapped **or < 50 % aligned** | a read spanning a sequence absent from the reference is clipped, not unmapped |
| short-read assembler | MEGAHIT | fast, good on uneven coverage |
| long-read assembler | Flye `--meta --min-overlap 1000` | uneven coverage of a read bin; reports circular contigs; Flye's automatic minimal overlap (≈ read N90) is too strict for read bins |
| HiFi tandem circles | on | small circles shorter than the reads are often missed by assemblers but appear as tandem copies inside reads |
| min contig length | 1 kb | shortest contigs in the training data |
| new contig | < 50 % aligned at ≥ 90 % identity (minimap2 `asm20`) | keeps divergent / chimeric contigs, removes genome copies |
| redundancy | contained ≥ 95 % at ≥ 95 % identity | Logan runs of one species overlap |
| Hi-C mapping | hicstuff, bowtie2, `cutsite` mode, MAPQ ≥ 30, PCR duplicates removed, **no** event filter | `--filter` removes the short-range pairs, i.e. all intra-contig contacts of small contigs (used for QC) |
| enzymes | from the protocol: Arima → `DpnII,HinfI`, Omni-C/Dovetail → `dnase`, Phase → `DpnII` | DToL Hi-C is mostly Arima v2 |
| profile resolution | 5 kb | resolution of the training data |
| min trans contacts | 500 | below, profiles are indistinguishable from Poisson noise |

---

## The classifier

The built-in model (`classification/models/hitchhiker_mlp.json`, numpy only) is the mean of five ReLU
networks (12 → 32 → 1) on **12 multi-scale features** of the contact profile against the host chromosomes:
log depth; excess variance and highest peak after aggregation in 5, 25, 100 and 500 kb windows; mean
obs/exp of the 1 % best 25-kb windows; share of contacts in the best 25-kb window; enrichment at
chromosome ends. Profiles are compared with the coverage-based expectation (bins with < 10 % of the median
coverage masked).

Training: positives = 2-micron (*S. cerevisiae*, 4 conditions), the 2-micron-like plasmid of *L.
fermentati*, Ddp5 (*D. discoideum*, 5 stages), with binomially thinned copies down to 500 contacts;
negatives = pRS416 (centromeric plasmid), pARS, mitochondria, rDNA, unplaced scaffolds, simulated
chromosomal fragments and Poisson null profiles. In leave-one-species-out validation (the predicted species
never seen in training) the AUC is 0.995; recovered contigs of 2-micron (*S. cerevisiae*) and Ddp5
(*D. discoideum*) are called STRONG_CANDIDATE with no false positive in those datasets.

**Limits.** Only three plasmid families were available for training: the model recognises plasmids that
hitchhike like 2-micron / Ddp5 and may miss molecules segregating differently. `chrom_end_oe` (enrichment
at chromosome ends) is the most informative feature and may partly reflect the organisation of the training
species. Treat STRONG_CANDIDATE as "worth a look" and check the plot and the QC fields. Another model with
the same JSON structure can be given with `--model`.

---

## Tests

```bash
pytest tests/test_units.py                 # parsers, DToL table, Logan headers, model (seconds)
bash tests/run_e2e.sh test_results         # end-to-end on simulated data (minutes, needs the conda env):
#   short reads + MEGAHIT, HiFi + Flye + tandem circles, mock Logan/ENA/SRA server, fake SLURM
python tests/real_classify.py out NAME file.mcool [NAME file.mcool ...]   # classify existing .mcool files
```

The simulated dataset (`tests/sim/make_sim.py`) contains four yeast chromosomes, a 6.3 kb circular plasmid,
a 30 kb piece removed from the reference and a contaminant bacterial genome. Checks: the plasmid is
recovered and circular (MEGAHIT end overlap, Flye, tandem copies, Logan self-loop), its relative copy number
is recovered from Logan, the missing piece is called CHROMOSOMAL_FRAGMENT at its true locus, and the
bacterial contigs are flagged as having no contact with the host. The simulated plasmid is a random
sequence with artificial, very sharp contact hotspots: it is **not expected** to be called STRONG_CANDIDATE
(it is not 2-micron-like); the classifier itself is checked on real contact maps with `tests/real_classify.py`.

---

## Changes from 0.0.x

* new commands `dtol`, `filter`, `run`, `status`; `classify` implemented; `pipeline` replaced by `run`
  (`--from/--to` instead of `--start/--stop`);
* one sample table for every command (DToL headers, old tables accepted), units = (species, reference);
* Logan over HTTP (no AWS CLI), reads of runs covered by Logan are not downloaded;
* reads: wget + local fasterq-dump (no prefetch), ENA fallback;
* long reads supported (minimap2 + Flye, tandem circles); contig re-alignment with minimap2;
* SLURM: real job chaining with dependencies, per-step resources, setup lines, waiting;
* hicstuff pairs → cooler matrices (hicstuff's internal cooler export is not used).

## References

1. Darwin Tree of Life — https://www.darwintreeoflife.org/
2. Logan — Chikhi R. et al. (2024) *Logan: planetary-scale genome assembly surveys life's diversity*, bioRxiv 10.1101/2024.07.30.605881
3. hicstuff — Matthey-Doret C. et al. (2020), Zenodo 10.5281/zenodo.4066363
4. cooler — Abdennur N., Mirny L. A. (2020) *Bioinformatics* 36:311
5. Girard F. et al. (2025) *Parasitic plasmids are anchored to inactive regions of eukaryotic chromosomes through a nucleosome signal*, EMBO J 44:2134
6. minimap2 — Li H. (2018) *Bioinformatics* 34:3094; MEGAHIT — Li D. et al. (2015) *Bioinformatics* 31:1674; Flye — Kolmogorov M. et al. (2020) *Nat. Methods* 17:1103

## Citing Kontiguity

Kontiguity, by M. Delouis, *(unpublished work)*, available at https://github.com/Mae-4815162342/kontiguity

## License

CC BY-NC 4.0 — see LICENSE.
