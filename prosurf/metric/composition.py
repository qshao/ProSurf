import numpy as np
import biotite.structure as struc

from prosurf.surface.sasa import residue_sasa, surface_residue_ids
from prosurf.surface.charges import assign_charges

_CHARGED = {"ASP", "GLU", "LYS", "ARG"}


def composition_features(arr, cfg) -> dict:
    """Bulk charge-composition features for one structure, independent of the
    spatial mixing captured by Z. These are the controls Z must beat.

    Returns floats:
      charged_frac           (#D,E,K,R residues) / (#residues), whole protein
      surface_charged_frac   surface charges / surface residues
      net_charge_per_res     |n_pos - n_neg| / (#residues), surface charges
      surface_charge_density surface charges / total surface SASA (Å^-2)
    """
    res_ids, res_names = struc.get_residues(arr)
    n_res = len(res_ids)
    if n_res == 0:
        return {"charged_frac": 0.0, "surface_charged_frac": 0.0,
                "net_charge_per_res": 0.0, "surface_charge_density": 0.0}

    n_charged = sum(1 for n in res_names if n in _CHARGED)
    charged_frac = n_charged / n_res

    surf_ids = surface_residue_ids(arr, cfg.rsasa_threshold)
    n_surf = len(surf_ids)
    charges = assign_charges(arr, surf_ids, his_weight=cfg.his_weight)
    n_pos = sum(1 for c in charges if c.sign > 0)
    n_neg = sum(1 for c in charges if c.sign < 0)
    n_surf_charged = n_pos + n_neg

    surface_charged_frac = (n_surf_charged / n_surf) if n_surf > 0 else 0.0
    surface_charged_frac = min(surface_charged_frac, 1.0)
    net_charge_per_res = abs(n_pos - n_neg) / n_res

    _, sasa = residue_sasa(arr)
    total_area = float(np.nansum(sasa))
    surface_charge_density = (n_surf_charged / total_area) if total_area > 0 else 0.0

    return {
        "charged_frac": float(charged_frac),
        "surface_charged_frac": float(surface_charged_frac),
        "net_charge_per_res": float(net_charge_per_res),
        "surface_charge_density": float(surface_charge_density),
    }
