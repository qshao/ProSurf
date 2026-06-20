from collections import namedtuple
import numpy as np

ProteinScore = namedtuple("ProteinScore", ["uniprot", "z_frac", "z_max", "z_mean", "n_patches"])

def normalize_density(raw_densities):
    cap = np.percentile(raw_densities, 99) if len(raw_densities) else 1.0
    cap = cap if cap > 0 else 1.0
    return lambda d: float(min(d / cap, 1.0))

def aggregate_protein(uniprot, patches, total_surface_area, n_surface_locations):
    if not patches:
        return ProteinScore(uniprot, 0.0, 0.0, 0.0, 0)
    sizes = np.array([p.size for p in patches])
    maxz = np.array([p.max_z for p in patches])
    meanz = np.array([p.mean_z for p in patches])
    z_frac = float(sizes.sum() / max(n_surface_locations, 1))
    z_max = float(maxz.max())
    z_mean = float(np.average(meanz, weights=sizes))
    return ProteinScore(uniprot, z_frac, z_max, z_mean, len(patches))
