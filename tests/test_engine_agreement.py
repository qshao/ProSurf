from pathlib import Path
from prosurf.config import MetricConfig
from prosurf.validate.robustness import engine_agreement


def test_engine_agreement_runs():
    cfg = MetricConfig(z_percentile=0.0)
    structs = [(Path("tests/data/mini.pdb"), "MINI"),
               (Path("tests/data/mini.pdb"), "MINI2")]
    rho = engine_agreement(structs, cfg)
    # identical inputs -> defined float (nan acceptable for degenerate set)
    assert isinstance(rho, float)
