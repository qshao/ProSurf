"""Tests for metric components: balance, density, mixing."""
import numpy as np
from prosurf.surface.charges import Charge
from prosurf.metric.components import balance, density, mixing_euclidean


def test_balance_perfect():
    assert balance(5, 5) == 1.0


def test_balance_one_sided():
    assert balance(5, 0) == 0.0


def test_balance_empty():
    assert balance(0, 0) == 0.0


def test_density():
    assert density(10.0, 100.0) == 0.1


def test_mixing_alternating_high():
    # + and - interleaved within 4 A -> high mixing
    charges = [
        Charge(1, 1, 1.0, np.array([0., 0., 0.])),
        Charge(2, -1, 1.0, np.array([3., 0., 0.])),
        Charge(3, 1, 1.0, np.array([6., 0., 0.])),
        Charge(4, -1, 1.0, np.array([9., 0., 0.])),
    ]
    assert mixing_euclidean(charges, neighbor_radius=4.0) > 0.9


def test_mixing_segregated_low():
    # all + clustered, all - far -> low mixing
    charges = [
        Charge(1, 1, 1.0, np.array([0., 0., 0.])),
        Charge(2, 1, 1.0, np.array([2., 0., 0.])),
        Charge(3, -1, 1.0, np.array([50., 0., 0.])),
        Charge(4, -1, 1.0, np.array([52., 0., 0.])),
    ]
    assert mixing_euclidean(charges, neighbor_radius=4.0) < 0.1
