import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra

_VDW = {"C": 1.7, "N": 1.55, "O": 1.52, "S": 1.8, "H": 1.2}


def _radii(arr):
    return np.array([_VDW.get(e, 1.7) for e in arr.element])


def sample_surface(arr, probe=1.4, density=1.0):
    coords = arr.coord
    radii = _radii(arr) + probe
    n_sphere = max(int(24 * density), 6)
    # Fibonacci sphere unit vectors
    i = np.arange(n_sphere)
    phi = np.arccos(1 - 2 * (i + 0.5) / n_sphere)
    theta = np.pi * (1 + 5**0.5) * i
    unit = np.column_stack([np.cos(theta) * np.sin(phi),
                            np.sin(theta) * np.sin(phi), np.cos(phi)])
    tree = cKDTree(coords)
    pts = []
    for c, r in zip(coords, radii):
        cand = c + unit * r
        # keep candidates not inside any other atom's probe sphere
        for p in cand:
            nbr = tree.query_ball_point(p, r)  # within own radius range
            inside = False
            for j in nbr:
                if np.linalg.norm(p - coords[j]) < radii[j] - 1e-6:
                    inside = True
                    break
            if not inside:
                pts.append(p)
    return np.array(pts) if pts else np.empty((0, 3))


def geodesic_graph(points, knn=8):
    tree = cKDTree(points)
    n = len(points)
    rows, cols, data = [], [], []
    k = min(knn + 1, n)
    dists, idx = tree.query(points, k=k)
    for i in range(n):
        for d, j in zip(dists[i][1:], idx[i][1:]):
            rows.append(i)
            cols.append(j)
            data.append(d)
    return csr_matrix((data, (rows, cols)), shape=(n, n))


def geodesic_distances(graph, source_idx):
    return dijkstra(graph, directed=False, indices=source_idx)
