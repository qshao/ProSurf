import numpy as np
from pathlib import Path
from prosurf.io.parse import load_structure
from prosurf.surface.sasa import residue_sasa, relative_sasa, surface_residue_ids

def test_residue_sasa_shape():
    arr = load_structure(Path("tests/data/mini.pdb"))
    res_ids, sasa = residue_sasa(arr)
    assert len(res_ids) == len(sasa) == 3
    assert np.all(sasa >= 0)

def test_relative_sasa_bounded_and_exposed():
    # isolated residues with no neighbors are highly exposed -> rsasa near/above many thresholds
    arr = load_structure(Path("tests/data/mini.pdb"))
    res_ids, rsasa = relative_sasa(arr)
    assert np.all(rsasa >= 0)
    assert np.all(rsasa[~np.isnan(rsasa)] > 0.20)  # all exposed in tiny structure

def test_surface_residue_ids():
    arr = load_structure(Path("tests/data/mini.pdb"))
    ids = surface_residue_ids(arr, threshold=0.20)
    assert len(ids) == 3
