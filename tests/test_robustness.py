import pandas as pd
from prosurf.validate.robustness import ranking_stability


def test_ranking_stability_identical_rankings():
    df = pd.DataFrame({
        "param_value": [0.1, 0.1, 0.2, 0.2],
        "uniprot": ["A", "B", "A", "B"],
        "z_mean": [0.9, 0.1, 0.8, 0.2],
    })
    assert abs(ranking_stability(df, "z_mean") - 1.0) < 1e-9


def test_ranking_stability_inverted():
    df = pd.DataFrame({
        "param_value": [0.1, 0.1, 0.2, 0.2],
        "uniprot": ["A", "B", "A", "B"],
        "z_mean": [0.9, 0.1, 0.1, 0.9],
    })
    assert ranking_stability(df, "z_mean") < 0
