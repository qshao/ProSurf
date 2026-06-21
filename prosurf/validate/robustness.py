from dataclasses import replace
from itertools import combinations

import pandas as pd
from scipy.stats import spearmanr

from prosurf.pipeline import analyze_structure


def sweep_parameter(structures, base_cfg, param, values):
    """Re-score a fixed set of proteins varying one MetricConfig field.

    Parameters
    ----------
    structures : list of (path, uniprot) tuples
    base_cfg   : MetricConfig
    param      : str — name of the MetricConfig field to vary
    values     : list — values to sweep over

    Returns
    -------
    pd.DataFrame with columns: uniprot, param_value, z_frac, z_max, z_mean, n_patches
    """
    rows = []
    for v in values:
        cfg = replace(base_cfg, **{param: v})
        for path, uniprot in structures:
            _, ps = analyze_structure(path, uniprot, cfg)
            row = ps._asdict()
            row["param_value"] = v
            rows.append(row)
    return pd.DataFrame(rows)


def ranking_stability(df, score_col="z_mean"):
    """Mean pairwise Spearman ρ of per-protein rankings across parameter values.

    Parameters
    ----------
    df       : pd.DataFrame with columns 'uniprot', 'param_value', and score_col
    score_col: str — column used as the score to rank proteins

    Returns
    -------
    float — mean Spearman ρ; returns 1.0 when there is only one param value
    """
    pivot = df.pivot(index="uniprot", columns="param_value", values=score_col)
    cols = list(pivot.columns)
    rhos = []
    for a, b in combinations(cols, 2):
        rho, _ = spearmanr(pivot[a], pivot[b])
        rhos.append(rho)
    return float(sum(rhos) / len(rhos)) if rhos else 1.0
