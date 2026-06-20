# prosurf/metric/engine_a.py
from collections import namedtuple
import numpy as np
from scipy.spatial import cKDTree
from prosurf.metric.components import balance, density, mixing_euclidean
from prosurf.surface.charges import assign_charges
from prosurf.surface.sasa import surface_residue_ids

LocationScore = namedtuple("LocationScore", ["res_id", "z", "n_pos", "n_neg", "xyz"])

def score_charges_a(charges, cfg):
    if not charges:
        return []
    coords = np.array([c.xyz for c in charges])
    signs = np.array([c.sign for c in charges])
    weights = np.array([c.weight for c in charges])
    tree = cKDTree(coords)
    disk_area = np.pi * cfg.patch_radius ** 2
    scores = []
    for i, c in enumerate(charges):
        idx = tree.query_ball_point(coords[i], cfg.patch_radius)
        region = [charges[j] for j in idx]
        n_pos = float(np.sum(weights[idx] * (signs[idx] > 0)))
        n_neg = float(np.sum(weights[idx] * (signs[idx] < 0)))
        B = balance(n_pos, n_neg)
        D = density(n_pos + n_neg, disk_area)
        M = mixing_euclidean(region, neighbor_radius=cfg.patch_radius / 2)
        scores.append(LocationScore(c.res_id, B * D * M, n_pos, n_neg, c.xyz))
    return scores  # density normalization across dataset applied in aggregation

def score_locations_a(arr, cfg):
    sids = surface_residue_ids(arr, cfg.rsasa_threshold)
    charges = assign_charges(arr, sids, his_weight=cfg.his_weight)
    return score_charges_a(charges, cfg)
