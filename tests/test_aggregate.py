import numpy as np
from prosurf.patches.cluster import Patch
from prosurf.patches.aggregate import aggregate_protein, normalize_density, ProteinScore

def test_normalize_density_caps_at_one():
    f = normalize_density(np.array([1., 2., 3., 100.]))
    assert 0.0 <= f(2.0) <= 1.0
    assert f(1000.0) == 1.0

def test_aggregate_protein_basic():
    patches = [
        Patch([1, 2], mean_z=0.5, max_z=0.6, n_pos=2, n_neg=2, size=2),
        Patch([5], mean_z=0.3, max_z=0.3, n_pos=1, n_neg=1, size=1),
    ]
    ps = aggregate_protein("P12345", patches, total_surface_area=100.0,
                           n_surface_locations=10)
    assert ps.uniprot == "P12345"
    assert ps.n_patches == 2
    assert ps.z_max == 0.6
    assert 0.0 <= ps.z_frac <= 1.0
    assert ps.z_mean > 0
