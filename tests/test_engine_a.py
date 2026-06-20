import numpy as np
from prosurf.config import MetricConfig
from prosurf.metric.engine_a import score_locations_a, LocationScore
from prosurf.surface.charges import Charge
import prosurf.metric.engine_a as ea

def test_score_locations_zwitterionic_beats_segregated(monkeypatch):
    cfg = MetricConfig(patch_radius=8.0)
    zwit = [Charge(i, (1 if i % 2 == 0 else -1), 1.0, np.array([float(i)*3, 0., 0.]))
            for i in range(6)]
    seg = [Charge(i, 1, 1.0, np.array([float(i)*3, 0., 0.])) for i in range(3)] + \
          [Charge(i, -1, 1.0, np.array([float(i)*3+50, 0., 0.])) for i in range(3)]
    z_zwit = max(s.z for s in ea.score_charges_a(zwit, cfg))
    z_seg = max(s.z for s in ea.score_charges_a(seg, cfg))
    assert z_zwit > z_seg

def test_score_charges_empty():
    cfg = MetricConfig()
    assert ea.score_charges_a([], cfg) == []
