import numpy as np
from scipy.spatial import cKDTree


def cross_colocation(coords, signs, d_max, geodesic_dist=None):
    """Fraction of short-range neighbor pairs that are opposite-signed,
    normalized against the expectation under random labeling.
    Returns ~1 when + and - co-locate (zwitterionic), ~0 when segregated."""
    n = len(coords)
    if n < 2:
        return 0.0
    pos_frac = np.mean(signs > 0)
    neg_frac = np.mean(signs < 0)
    expected_opposite = 2 * pos_frac * neg_frac  # random-labeling baseline
    if expected_opposite == 0:
        return 0.0
    if geodesic_dist is None:
        tree = cKDTree(coords)
        pairs = tree.query_pairs(d_max, output_type="ndarray")
    else:
        pairs = np.array([[i, j] for i in range(n) for j in range(i+1, n)
                          if geodesic_dist(i, j) <= d_max])
    if len(pairs) == 0:
        return 0.0
    opp = np.mean(signs[pairs[:, 0]] != signs[pairs[:, 1]])
    return float(min(opp / expected_opposite, 1.0))
