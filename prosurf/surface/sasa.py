import numpy as np
import biotite.structure as struc

# Tien et al. 2013 theoretical max ASA (Angstrom^2)
MAX_ASA = {
    "ALA":129.0,"ARG":274.0,"ASN":195.0,"ASP":193.0,"CYS":167.0,"GLU":223.0,
    "GLN":225.0,"GLY":104.0,"HIS":224.0,"ILE":197.0,"LEU":201.0,"LYS":236.0,
    "MET":224.0,"PHE":240.0,"PRO":159.0,"SER":155.0,"THR":172.0,"TRP":285.0,
    "TYR":263.0,"VAL":174.0,
}

def residue_sasa(arr):
    atom_sasa = struc.sasa(arr, vdw_radii="ProtOr")
    res_ids, _ = struc.get_residues(arr)
    sasa = struc.apply_residue_wise(arr, atom_sasa, np.nansum)
    return res_ids, sasa

def relative_sasa(arr):
    res_ids, sasa = residue_sasa(arr)
    _, res_names = struc.get_residues(arr)
    maxasa = np.array([MAX_ASA.get(n, np.nan) for n in res_names])
    return res_ids, sasa / maxasa

def surface_residue_ids(arr, threshold):
    res_ids, rsasa = relative_sasa(arr)
    return res_ids[np.nan_to_num(rsasa, nan=0.0) > threshold]
