from pathlib import Path
from prosurf.config import MetricConfig
from prosurf.io.parse import load_structure
from prosurf.metric.composition import composition_features


def test_composition_on_mini():
    # mini.pdb = LYS, ASP, ALA -> 2 of 3 residues charged
    arr = load_structure(Path("tests/data/mini.pdb"))
    cfg = MetricConfig()
    feats = composition_features(arr, cfg)
    assert abs(feats["charged_frac"] - 2 / 3) < 1e-6
    for k in ("surface_charged_frac", "net_charge_per_res", "surface_charge_density"):
        assert feats[k] >= 0.0
    assert 0.0 <= feats["surface_charged_frac"] <= 1.0


def test_composition_keys_are_floats():
    arr = load_structure(Path("tests/data/mini.pdb"))
    feats = composition_features(arr, MetricConfig())
    assert set(feats) == {"charged_frac", "surface_charged_frac",
                          "net_charge_per_res", "surface_charge_density"}
    assert all(isinstance(v, float) for v in feats.values())
