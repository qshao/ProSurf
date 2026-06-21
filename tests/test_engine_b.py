# tests/test_engine_b.py
from pathlib import Path
from prosurf.config import MetricConfig
from prosurf.pipeline import analyze_structure


def test_engine_b_runs_and_shares_shape():
    cfg = MetricConfig(z_percentile=0.0)
    patches_a, score_a = analyze_structure(Path("tests/data/mini.pdb"), "MINI", cfg, engine="a")
    patches_b, score_b = analyze_structure(Path("tests/data/mini.pdb"), "MINI", cfg, engine="b")
    assert score_a._fields == score_b._fields
    assert score_b.uniprot == "MINI"
