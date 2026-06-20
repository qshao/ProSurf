"""Metric components: balance, density, and mixing."""
import numpy as np
from scipy.spatial import cKDTree


def balance(n_pos, n_neg):
    """Compute neutrality score (balance) between 0 and 1.

    Args:
        n_pos: Number or count of positive charges.
        n_neg: Number or count of negative charges.

    Returns:
        float: Neutrality score. 1.0 when perfectly balanced (equal positive/negative),
               0.0 when completely one-sided or empty.
    """
    total = n_pos + n_neg
    if total == 0:
        return 0.0
    return 1.0 - abs(n_pos - n_neg) / total


def density(total_charge, area):
    """Compute raw charge density.

    This is the raw, un-normalized density. Normalization across the dataset
    happens at a later stage (e.g., in engine or aggregation).

    Args:
        total_charge: Sum of all charge weights on the surface.
        area: Surface area (e.g., SASA in Angstrom^2).

    Returns:
        float: Charge density (total_charge / area). Returns 0.0 if area <= 0.
    """
    if area <= 0:
        return 0.0
    return total_charge / area


def mixing_euclidean(charges, neighbor_radius):
    """Compute fraction of opposite-sign neighbors within Euclidean distance.

    For each charge, finds all neighbors within neighbor_radius (Euclidean distance),
    then computes the fraction of those neighbors with opposite sign. The overall
    mixing is the mean fraction across all charges.

    Args:
        charges: list[Charge] — namedtuples with .xyz (coordinate) and .sign fields.
        neighbor_radius: float — radius (Angstroms) for neighbor search.

    Returns:
        float: Mean fraction of opposite-sign neighbors (0.0 to 1.0).
               Returns 0.0 if fewer than 2 charges or no neighbors found.
    """
    if len(charges) < 2:
        return 0.0

    coords = np.array([c.xyz for c in charges])
    signs = np.array([c.sign for c in charges])

    tree = cKDTree(coords)
    fracs = []

    for i, c in enumerate(charges):
        idx = tree.query_ball_point(coords[i], neighbor_radius)
        idx = [j for j in idx if j != i]  # Exclude self

        if not idx:
            continue

        opposite = np.sum(signs[idx] != signs[i])
        fracs.append(opposite / len(idx))

    return float(np.mean(fracs)) if fracs else 0.0
