# prosurf/metric/engine_a.py
from collections import namedtuple
import numpy as np
from scipy.spatial import cKDTree
from prosurf.metric.components import balance, density, mixing_euclidean
from prosurf.surface.charges import assign_charges
from prosurf.surface.sasa import surface_residue_ids

LocationScore = namedtuple("LocationScore", ["res_id", "z", "n_pos", "n_neg", "xyz"])

_MIN_PATCH_RADIUS = 8.0  # Å floor so mixing sees at least a few neighbors

def score_charges_a(charges, cfg):
    if not charges:
        return []
    coords = np.array([c.xyz for c in charges])
    signs = np.array([c.sign for c in charges])
    weights = np.array([c.weight for c in charges])
    tree = cKDTree(coords)

    # Adaptive radius: scale with the spatial extent of the charge cloud
    centroid = coords.mean(axis=0)
    R_g = float(np.sqrt(((coords - centroid) ** 2).sum(axis=1).mean()))
    r = max(cfg.patch_radius_frac * R_g, _MIN_PATCH_RADIUS)
    disk_area = np.pi * r ** 2

    # First pass: collect (n_pos, n_neg, region, idx) per location
    per_location = []
    for i, c in enumerate(charges):
        idx = tree.query_ball_point(coords[i], r)
        region = [charges[j] for j in idx]
        n_pos = float(np.sum(weights[idx] * (signs[idx] > 0)))
        n_neg = float(np.sum(weights[idx] * (signs[idx] < 0)))
        per_location.append((c, idx, region, n_pos, n_neg))

    # Compute raw densities and normalize (D̂ ∈ [0, 1])
    raw_densities = np.array([density(n_pos + n_neg, disk_area)
                               for _, _, _, n_pos, n_neg in per_location])
    d_max = raw_densities.max() if raw_densities.max() > 0 else 1.0

    scores = []
    for (c, idx, region, n_pos, n_neg), raw_d in zip(per_location, raw_densities):
        B = balance(n_pos, n_neg)
        D_hat = raw_d / d_max
        M = mixing_euclidean(region, neighbor_radius=r / 2)
        scores.append(LocationScore(c.res_id, B * D_hat * M, n_pos, n_neg, c.xyz))
    return scores

def score_locations_a(arr, cfg):
    sids = surface_residue_ids(arr, cfg.rsasa_threshold)
    charges = assign_charges(arr, sids, his_weight=cfg.his_weight)
    return score_charges_a(charges, cfg)
