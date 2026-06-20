# tests/test_cluster.py
import numpy as np
from prosurf.config import MetricConfig
from prosurf.metric.engine_a import LocationScore
from prosurf.patches.cluster import cluster_patches


def test_two_separate_patches():
    cfg = MetricConfig(z_percentile=0.0)  # keep all
    scores = [
        LocationScore(1, 0.9, 2, 2, np.array([0., 0., 0.])),
        LocationScore(2, 0.8, 2, 2, np.array([3., 0., 0.])),
        LocationScore(3, 0.9, 2, 2, np.array([100., 0., 0.])),
    ]
    patches = cluster_patches(scores, cfg, adjacency_radius=8.0)
    assert len(patches) == 2
    sizes = sorted(p.size for p in patches)
    assert sizes == [1, 2]


def test_percentile_filters_low_z():
    cfg = MetricConfig(z_percentile=90.0)
    scores = [LocationScore(i, float(i)/10, 1, 1, np.array([float(i), 0., 0.]))
              for i in range(11)]
    patches = cluster_patches(scores, cfg, adjacency_radius=8.0)
    kept = sum(p.size for p in patches)
    assert kept <= 2  # only top ~10%
