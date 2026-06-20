import numpy as np
from pathlib import Path
from prosurf.io.parse import load_structure
from prosurf.surface.charges import assign_charges, Charge

def test_assign_charges_signs():
    arr = load_structure(Path("tests/data/mini.pdb"))
    surface_ids = np.array([1, 2, 3])
    charges = assign_charges(arr, surface_ids, his_weight=0.0)
    signs = {c.res_id: c.sign for c in charges if c.weight == 1.0}
    # Lys is +, Asp is -, Ala has no side-chain charge
    lys_id = [c.res_id for c in charges if c.sign == 1 and c.weight == 1.0]
    asp_id = [c.res_id for c in charges if c.sign == -1 and c.weight == 1.0]
    assert len(lys_id) >= 1
    assert len(asp_id) >= 1

def test_assign_charges_excludes_buried():
    arr = load_structure(Path("tests/data/mini.pdb"))
    charges = assign_charges(arr, surface_ids=np.array([], dtype=int), his_weight=0.0)
    # no surface residues -> only termini charges remain
    sidechain = [c for c in charges if c.weight == 1.0 and c.sign != 0]
    assert all(c.res_id in (1, 3) for c in sidechain) or len(sidechain) == 0
