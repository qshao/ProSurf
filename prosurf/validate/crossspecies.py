import numpy as np
from scipy import stats


def partial_spearman(x, y, covariates):
    """Partial Spearman correlation of x and y controlling for covariates.

    Rank-transforms all variables, regresses ranked x and ranked y on the
    ranked covariates (with intercept), and correlates the residuals.
    `covariates` is a list of 1-D array-likes. Returns (partial_r, p_value).
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    C = np.column_stack([np.asarray(c, float) for c in covariates])

    xr = stats.rankdata(x)
    yr = stats.rankdata(y)
    Cr = np.column_stack([stats.rankdata(C[:, j]) for j in range(C.shape[1])])
    A = np.column_stack([np.ones(len(xr)), Cr])

    bx, *_ = np.linalg.lstsq(A, xr, rcond=None)
    by, *_ = np.linalg.lstsq(A, yr, rcond=None)
    rx = xr - A @ bx
    ry = yr - A @ by
    r, p = stats.pearsonr(rx, ry)
    return float(r), float(p)


def nested_regression_delta(df, target, base_features, added_feature):
    """Compare OLS models target ~ base vs target ~ base + added_feature.

    Predictors are z-standardized. Returns dict: r2_base, r2_full, delta_r2,
    f_stat, p_value (F-test for the single added term), added_coef (standardized).
    """
    y = df[target].to_numpy(float)
    n = len(y)

    def standardize(M):
        M = np.asarray(M, float)
        return (M - M.mean(0)) / M.std(0)

    Xb = standardize(df[base_features].to_numpy(float))
    Xf = standardize(df[base_features + [added_feature]].to_numpy(float))
    Ab = np.column_stack([np.ones(n), Xb])
    Af = np.column_stack([np.ones(n), Xf])

    def fit(A):
        beta, *_ = np.linalg.lstsq(A, y, rcond=None)
        resid = y - A @ beta
        rss = float(resid @ resid)
        tss = float(((y - y.mean()) ** 2).sum())
        return beta, rss, 1.0 - rss / tss

    bb, rss_b, r2_b = fit(Ab)
    bf, rss_f, r2_f = fit(Af)

    df_num = Af.shape[1] - Ab.shape[1]
    df_den = n - Af.shape[1]
    f_stat = ((rss_b - rss_f) / df_num) / (rss_f / df_den)
    p_value = float(stats.f.sf(f_stat, df_num, df_den))

    return {
        "r2_base": r2_b, "r2_full": r2_f, "delta_r2": r2_f - r2_b,
        "f_stat": float(f_stat), "p_value": p_value, "added_coef": float(bf[-1]),
    }
