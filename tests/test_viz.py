from pathlib import Path
import matplotlib
matplotlib.use("Agg")
from prosurf.config import MetricConfig
from prosurf.io.parse import load_structure
from prosurf.metric.engine_a import score_locations_a
from prosurf.viz.patchmap import render_patch_map


def test_render_creates_png(tmp_path):
    cfg = MetricConfig(z_percentile=0.0)
    arr = load_structure(Path("tests/data/mini.pdb"))
    scores = score_locations_a(arr, cfg)
    out = render_patch_map(arr, scores, tmp_path / "map.png")
    assert out.exists() and out.stat().st_size > 0
