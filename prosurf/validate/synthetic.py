import numpy as np
from prosurf.surface.charges import Charge

def make_synthetic(pattern, n, spacing=3.0):
    charges = []
    for i in range(n):
        xyz = np.array([i * spacing, 0.0, 0.0])
        if pattern == "alternating":
            sign = 1 if i % 2 == 0 else -1
        elif pattern == "segregated":
            sign = 1 if i < n // 2 else -1
        elif pattern == "charge_free":
            continue
        else:
            raise ValueError(pattern)
        charges.append(Charge(i, sign, 1.0, xyz))
    return charges
