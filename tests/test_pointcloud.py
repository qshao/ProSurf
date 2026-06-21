import numpy as np
from pathlib import Path
from prosurf.io.parse import load_structure
from prosurf.surface.pointcloud import sample_surface, geodesic_graph, geodesic_distances

def test_sample_surface_nonempty():
    arr = load_structure(Path("tests/data/mini.pdb"))
    pts = sample_surface(arr, density=0.5)
    assert pts.ndim == 2 and pts.shape[1] == 3 and len(pts) > 0

def test_geodesic_distances_monotone():
    pts = np.array([[0,0,0],[1,0,0],[2,0,0],[3,0,0]], dtype=float)
    g = geodesic_graph(pts, knn=2)
    d = geodesic_distances(g, 0)
    assert d[0] == 0 and d[3] >= d[1]
