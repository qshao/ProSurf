import numpy as np
from prosurf.metric.mixing_spatial import cross_colocation


def test_colocation_alternating_high():
    coords = np.array([[i*3.,0,0] for i in range(8)])
    signs = np.array([1,-1,1,-1,1,-1,1,-1])
    assert cross_colocation(coords, signs, d_max=4.0) > 0.8


def test_colocation_segregated_low():
    coords = np.array([[i*3.,0,0] for i in range(8)])
    signs = np.array([1,1,1,1,-1,-1,-1,-1])
    assert cross_colocation(coords, signs, d_max=4.0) < 0.4
