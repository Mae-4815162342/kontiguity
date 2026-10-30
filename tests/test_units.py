"""Fast unit tests (no external tool, no network):  pytest tests/test_units.py"""
import json
import os

import numpy as np
import pandas as pd
import pytest

from kontiguity import dtol, logan, samples
from kontiguity.fasta import parse_ena_header
from kontiguity.classification.model import MLPEnsemble
from kontiguity.classification.features import hand_feats, NAMES

HERE = os.path.dirname(__file__)
DATA = os.path.join(HERE, "..", "test_data", "load")


def test_ena_header():
    assert parse_ena_header("ENA|OX637620|OX637620.1 Tausonia pullulans genome assembly, chromosome: 1")[:3] == \
        ("OX637620.1", "chromosome", "1")
    assert parse_ena_header("ENA|OX6|OX6.1 Some species mitochondrion, complete genome")[1] == "organelle"
    assert parse_ena_header("NC_001133.9 chrI")[0] == "NC_001133.9"


def test_logan_header():
    h = logan.parse_header("SRR1_12 ka:f:224.7 L:+:12:+ L:-:3:+")
    assert h["circular"] and h["ka"] == pytest.approx(224.7) and h["index"] == "12"
    assert not logan.parse_header("SRR1_3 ka:f:2 L:+:3:-")["circular"]   # hairpin, not a circle
    assert not logan.parse_header("SRR1_3 km:f:2 L:+:4:+")["circular"]


def test_inference():
    assert samples.infer_tech("PacBio - HiFi") == "hifi"
    assert samples.infer_tech("Chromium genome", layout="PAIRED") == "short"
    assert samples.infer_tech("", "OXFORD_NANOPORE PromethION") == "ont"
    assert samples.infer_enzymes("Hi-C - Arima v2") == "DpnII,HinfI"
    assert samples.infer_enzymes("Dovetail Omni-C") == "dnase"


def test_dtol_table_grouping(tmp_path):
    """The legacy one-row-per-combination DToL table is grouped into (species, reference) units."""
    df = samples.to_frame(samples.rows_from_table(os.path.join(DATA, "test_dataset.csv")))
    hm = df[df.name == "Hydraecia micacea"].iloc[0]
    assert hm.wgs.split(";") == ["ERR6688535", "ERR6688537", "ERR6688538", "ERR6688536", "ERR6939243"]
    assert hm.wgs_tech.split(";") == ["short"] * 4 + ["hifi"]
    assert hm.hic == "ERR6688534" and hm.enzymes == "DpnII,HinfI"
    vv = df[df.name == "Vespula vulgaris"].iloc[0]
    assert vv.hic == "ERR6054674;ERR6054675" and vv.unit == "Vespula_vulgaris/GCA_905475345.1"


def test_direct_inputs(tmp_path):
    r1, r2 = tmp_path / "x_R1.fq.gz", tmp_path / "x_R2.fq.gz"
    r1.write_bytes(b""); r2.write_bytes(b"")
    df = samples.to_frame(samples.rows_from_args("My species", "GCA_000001.1,GCA_000002.1", f"{r1}:{r2},SRR123", "ERR9",
                                                 wgs_tech="short,hifi"))
    assert list(df.unit) == ["My_species/GCA_000001.1", "My_species/GCA_000002.1"]
    rows = samples.reads_of(df.iloc[0], "wgs")
    assert rows[0]["layout"] == "PAIRED" and len(rows[0]["files"]) == 2 and rows[1]["id"] == "SRR123"
    assert rows[1]["tech"] == "hifi"


def test_dtol_records():
    rec = [{"_id": "Vespula vulgaris", "_source": {
        "tax_id": "7454",
        "experiment": [
            {"library_strategy": "WGS", "run_accession": "ERR1", "library_layout": "SINGLE",
             "library_construction_protocol": "PacBio - HiFi", "study_accession": "PRJEB1", "instrument_platform": "PACBIO_SMRT"},
            {"library_strategy": "Hi-C", "run_accession": "ERR2", "library_layout": "PAIRED",
             "library_construction_protocol": "Hi-C - Arima v2", "study_accession": "PRJEB1"},
            {"library_strategy": "RNA-Seq", "run_accession": "ERR3"}],
        "assemblies": [{"accession": "GCA_905475345", "version": "1", "assembly_name": "iyVesVulg1.1"},
                       {"accession": "GCA_905475355.1", "assembly_name": "iyVesVulg1.1 alternate haplotype"}]}},
        {"_id": "No hic species", "_source": {"experiment": [{"library_strategy": "WGS", "run_accession": "ERR9"}],
                                               "assemblies": [{"accession": "GCA_1.1"}]}}]
    t = dtol.records_to_table(rec)
    assert list(t.ref) == ["GCA_905475345.1", "GCA_905475355.1"] and t.wgs.iloc[0] == "ERR1"
    assert len(dtol.records_to_table(rec, primary_only=True)) == 1
    df = samples.to_frame(samples.rows_from_table(_save(t)))
    assert df.wgs_tech.iloc[0] == "hifi"


def _save(t):
    import tempfile
    p = tempfile.mktemp(suffix=".tsv")
    t.to_csv(p, sep="\t", index=False)
    return p


def test_model_and_features():
    m = MLPEnsemble()
    assert m.features == NAMES
    rng = np.random.default_rng(0)
    n = 2000
    w = np.ones(n) / n
    chrom = np.repeat(np.arange(10), n // 10)
    flat = rng.multinomial(20000, w)[None]
    f = hand_feats(flat, w, np.ones(n), chrom)
    assert f.shape == (1, 12)
    assert 0 <= m.predict_proba(f)[0] <= 1


def test_platform_used_when_protocol_is_generic(tmp_path):
    p = tmp_path / "t.tsv"
    p.write_text("name\tref\twgs\tlayout_WGS\tprotocol_WGS\tplatform_WGS\thic\n"
                 "A b\tGCA_1.1\tERR1;ERR2;ERR3\tSINGLE;SINGLE;PAIRED\tWGS;WGS;Chromium genome\t"
                 "OXFORD_NANOPORE PromethION;PACBIO_SMRT Sequel;ILLUMINA\tERR9\n")
    df = samples.to_frame(samples.rows_from_table(str(p)))
    assert df.wgs_tech.iloc[0] == "ont;clr;short"


def test_tmp_dir_layout(tmp_path, monkeypatch):
    """--tmp-dir moves the bulky files out of the output tree; '$VAR' is expanded at run time."""
    from kontiguity.utils.runner import Layout, set_tmp_dir
    out = tmp_path / "out"
    row = dict(species="Sp_a", unit="Sp_a/REF", ref_slug="REF", name="Sp a")
    set_tmp_dir(str(out), None)
    lay = Layout(str(out), row)
    assert lay.wgs_dir.startswith(str(out)) and lay.mcool_store == lay.mcool
    monkeypatch.setenv("MYSCRATCH", str(tmp_path / "scratch"))
    root = set_tmp_dir(str(out), "$MYSCRATCH/kt")
    lay = Layout(str(out), row)
    assert root.startswith(str(tmp_path / "scratch" / "kt" / "kontiguity_out_"))
    for p in (lay.wgs_dir, lay.hic_dir_reads, lay.logan_dir, lay.retrieve_work, lay.hic_work, lay.mcool_store):
        assert p.startswith(root), p
    for p in (lay.genome_fa, lay.report, lay.mcool, lay.raw_contigs):
        assert p.startswith(str(out)), p
    set_tmp_dir(str(out), None)
