# prosurf/patches/cluster.py
from collections import namedtuple
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

Patch = namedtuple("Patch", ["res_ids", "mean_z", "max_z", "n_pos", "n_neg", "size"])

def cluster_patches(scores, cfg, adjacency_radius=8.0):
    """
    Cluster locations into patches using connected components.

    Parameters
    ----------
    scores : list[LocationScore]
        List of LocationScore namedtuples with res_id, z, n_pos, n_neg, xyz.
    cfg : MetricConfig
        Configuration object with z_percentile field.
    adjacency_radius : float, optional
        Radius within which locations are considered adjacent (default: 8.0).

    Returns
    -------
    list[Patch]
        List of Patch namedtuples, one per connected component.
    """
    if not scores:
        return []

    # Filter locations: keep only those with z >= percentile AND z > 0
    z = np.array([s.z for s in scores])
    cutoff = np.percentile(z, cfg.z_percentile)
    keep = [i for i, s in enumerate(scores) if s.z >= cutoff and s.z > 0]

    if not keep:
        return []

    # Build adjacency graph from kept locations
    kept = [scores[i] for i in keep]
    coords = np.array([s.xyz for s in kept])
    tree = cKDTree(coords)
    pairs = tree.query_pairs(adjacency_radius, output_type="ndarray")

    # Create sparse adjacency matrix
    n = len(kept)
    if len(pairs) > 0:
        data = np.ones(len(pairs))
        graph = csr_matrix((data, (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    else:
        graph = csr_matrix((n, n))

    # Find connected components
    n_comp, labels = connected_components(graph, directed=False)

    # Create patches from components
    patches = []
    for c in range(n_comp):
        members = [kept[i] for i in range(n) if labels[i] == c]
        zs = np.array([m.z for m in members])
        patches.append(Patch(
            res_ids=[m.res_id for m in members],
            mean_z=float(zs.mean()),
            max_z=float(zs.max()),
            n_pos=float(sum(m.n_pos for m in members)),
            n_neg=float(sum(m.n_neg for m in members)),
            size=len(members)))

    return patches
