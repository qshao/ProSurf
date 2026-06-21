import numpy as np
import pandas as pd
from prosurf.validate.crossspecies import partial_spearman, nested_regression_delta


def test_partial_spearman_removes_confound():
    rng = np.random.default_rng(0)
    z = rng.normal(size=400)              # confound
    x = z + 0.01 * rng.normal(size=400)   # x driven only by z
    y = z + 0.01 * rng.normal(size=400)   # y driven only by z
    # raw correlation high, partial (controlling z) near zero
    r_partial, p = partial_spearman(x, y, [z])
    assert abs(r_partial) < 0.2


def test_partial_spearman_keeps_real_signal():
    rng = np.random.default_rng(1)
    z = rng.normal(size=400)
    sig = rng.normal(size=400)
    x = z + sig
    y = z + sig                            # shared signal beyond z
    r_partial, p = partial_spearman(x, y, [z])
    assert r_partial > 0.4
    assert p < 1e-3


def test_nested_regression_delta_detects_added_signal():
    rng = np.random.default_rng(2)
    base = rng.normal(size=500)
    extra = rng.normal(size=500)
    y = 1.0 * base + 1.0 * extra + 0.1 * rng.normal(size=500)
    df = pd.DataFrame({"y": y, "base": base, "extra": extra})
    res = nested_regression_delta(df, "y", ["base"], "extra")
    assert res["delta_r2"] > 0.2
    assert res["p_value"] < 1e-6
    assert res["added_coef"] > 0


def test_nested_regression_delta_null_added_signal():
    rng = np.random.default_rng(3)
    base = rng.normal(size=500)
    noise = rng.normal(size=500)           # unrelated to y
    y = base + 0.1 * rng.normal(size=500)
    df = pd.DataFrame({"y": y, "base": base, "noise": noise})
    res = nested_regression_delta(df, "y", ["base"], "noise")
    assert res["delta_r2"] < 0.02
    assert res["p_value"] > 0.05
