from collections import namedtuple
import numpy as np
import biotite.structure as struc

Charge = namedtuple("Charge", ["res_id", "sign", "weight", "xyz"])

def _atom_coord(arr, res_id, atom_name):
    mask = (arr.res_id == res_id) & (arr.atom_name == atom_name)
    return arr.coord[mask][0] if mask.any() else None

def _center(arr, res_id, names):
    coords = [_atom_coord(arr, res_id, n) for n in names]
    coords = [c for c in coords if c is not None]
    return np.mean(coords, axis=0) if coords else None

def assign_charges(arr, surface_ids, his_weight=0.0):
    surface = set(int(i) for i in surface_ids)
    charges = []
    res_ids, res_names = struc.get_residues(arr)
    for rid, rname in zip(res_ids, res_names):
        if int(rid) not in surface:
            continue
        if rname == "LYS":
            xyz = _atom_coord(arr, rid, "NZ")
            if xyz is not None: charges.append(Charge(int(rid), 1, 1.0, xyz))
        elif rname == "ARG":
            xyz = _atom_coord(arr, rid, "CZ")
            if xyz is not None: charges.append(Charge(int(rid), 1, 1.0, xyz))
        elif rname == "ASP":
            xyz = _center(arr, rid, ["OD1", "OD2"])
            if xyz is not None: charges.append(Charge(int(rid), -1, 1.0, xyz))
        elif rname == "GLU":
            xyz = _center(arr, rid, ["OE1", "OE2"])
            if xyz is not None: charges.append(Charge(int(rid), -1, 1.0, xyz))
        elif rname == "HIS" and his_weight > 0:
            xyz = _center(arr, rid, ["ND1", "NE2"])
            if xyz is not None: charges.append(Charge(int(rid), 1, his_weight, xyz))
    # termini as point charges at terminal CA-adjacent atoms
    # Skip terminus charge if a same-sign sidechain charge is already assigned
    # (e.g. N-terminus +1 is redundant if the first residue is LYS or ARG)
    nterm, cterm = res_ids[0], res_ids[-1]
    first_res_has_pos = any(c.res_id == int(nterm) and c.sign == 1 for c in charges)
    last_res_has_neg = any(c.res_id == int(cterm) and c.sign == -1 for c in charges)
    n_xyz = _atom_coord(arr, nterm, "N")
    c_xyz = _atom_coord(arr, cterm, "C")
    if n_xyz is not None and not first_res_has_pos:
        charges.append(Charge(int(nterm), 1, 1.0, n_xyz))
    if c_xyz is not None and not last_res_has_neg:
        charges.append(Charge(int(cterm), -1, 1.0, c_xyz))
    return charges
