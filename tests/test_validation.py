from prosurf.config import MetricConfig
from prosurf.validate.synthetic import make_synthetic
from prosurf.validate.controls import protein_zwit_score, auroc

def test_synthetic_ordering():
    cfg = MetricConfig(patch_radius=8.0)
    alt = protein_zwit_score(make_synthetic("alternating", 12), cfg)
    seg = protein_zwit_score(make_synthetic("segregated", 12), cfg)
    free = protein_zwit_score(make_synthetic("charge_free", 12), cfg)
    assert alt > seg
    assert alt > free

def test_auroc_perfect_separation():
    assert auroc([0.9, 0.8, 0.95], [0.1, 0.2, 0.05]) == 1.0

def test_auroc_random():
    assert abs(auroc([0.5, 0.5], [0.5, 0.5]) - 0.5) < 1e-9
