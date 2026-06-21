# prosurf/metric/engine_b.py
import numpy as np
from scipy.spatial import cKDTree
from prosurf.metric.engine_a import LocationScore
from prosurf.metric.components import balance, density
from prosurf.metric.mixing_spatial import cross_colocation
from prosurf.surface.charges import assign_charges
from prosurf.surface.sasa import surface_residue_ids


def score_locations_b(arr, cfg):
    sids = surface_residue_ids(arr, cfg.rsasa_threshold)
    charges = assign_charges(arr, sids, his_weight=cfg.his_weight)
    if not charges:
        return []
    coords = np.array([c.xyz for c in charges])
    signs = np.array([c.sign for c in charges])
    weights = np.array([c.weight for c in charges])
    tree = cKDTree(coords)
    disk_area = np.pi * cfg.patch_radius ** 2

    # First pass: collect per-location counts and indices
    per_location = []
    for i, c in enumerate(charges):
        idx = tree.query_ball_point(coords[i], cfg.patch_radius)
        n_pos = float(np.sum(weights[idx] * (signs[idx] > 0)))
        n_neg = float(np.sum(weights[idx] * (signs[idx] < 0)))
        per_location.append((c, idx, n_pos, n_neg))

    # Compute raw densities and normalize (D̂ ∈ [0, 1])
    raw_densities = np.array([density(n_pos + n_neg, disk_area)
                               for _, _, n_pos, n_neg in per_location])
    d_max = raw_densities.max() if raw_densities.max() > 0 else 1.0

    scores = []
    for (c, idx, n_pos, n_neg), raw_d in zip(per_location, raw_densities):
        B = balance(n_pos, n_neg)
        D_hat = raw_d / d_max
        M = cross_colocation(coords[idx], signs[idx], d_max=cfg.patch_radius / 2)
        scores.append(LocationScore(c.res_id, B * D_hat * M, n_pos, n_neg, c.xyz))
    return scores
