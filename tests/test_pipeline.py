# tests/test_pipeline.py
from pathlib import Path
from prosurf.config import MetricConfig
from prosurf.pipeline import analyze_structure

def test_analyze_structure_runs_on_mini():
    cfg = MetricConfig(z_percentile=0.0)
    patches, score = analyze_structure(Path("tests/data/mini.pdb"), "MINI", cfg)
    assert score.uniprot == "MINI"
    assert isinstance(patches, list)
    assert score.n_patches >= 0
