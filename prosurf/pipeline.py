# prosurf/pipeline.py
import pickle
import numpy as np
import pandas as pd
from pathlib import Path
from prosurf.io.parse import load_structure
from prosurf.io.fetch import fetch_af2
from prosurf.metric.engine_a import score_locations_a
from prosurf.metric.engine_b import score_locations_b
from prosurf.patches.cluster import cluster_patches
from prosurf.patches.aggregate import aggregate_protein
from prosurf.surface.sasa import residue_sasa, surface_residue_ids


def analyze_structure(path, uniprot, cfg, engine="a"):
    """
    Analyze a single structure and return patches and protein-level score.

    Parameters
    ----------
    path : Path
        Path to the PDB/CIF file.
    uniprot : str
        UniProt accession (used as identifier in the score).
    cfg : MetricConfig
        Metric configuration (thresholds, weights, etc.).
    engine : str, optional
        Scoring engine to use. Currently only "a" is supported (default: "a").
        Task 16 will add engine "b" and a selector here.

    Returns
    -------
    tuple[list[Patch], ProteinScore]
    """
    arr = load_structure(path)

    scorer = score_locations_a if engine == "a" else score_locations_b
    scores = scorer(arr, cfg)

    patches = cluster_patches(scores, cfg)
    _, sasa = residue_sasa(arr)
    total_area = float(np.nansum(sasa))
    n_locations = len(surface_residue_ids(arr, cfg.rsasa_threshold))
    ps = aggregate_protein(uniprot, patches, total_area, max(n_locations, 1))
    return patches, ps


def run_pilot(uniprots, paths, cfg):
    """
    Fetch AF2 structures for each UniProt ID, analyze them, cache patch lists,
    and return a DataFrame with one row per protein.

    Parameters
    ----------
    uniprots : list[str]
        UniProt accession IDs to process.
    paths : PathsConfig
        Paths for data and cache directories.
    cfg : MetricConfig
        Metric configuration.

    Returns
    -------
    pandas.DataFrame
        One row per protein with ProteinScore fields as columns.
    """
    paths.cache_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for up in uniprots:
        pdb = fetch_af2(up, paths.data_dir)
        patches, ps = analyze_structure(pdb, up, cfg)
        with open(paths.cache_dir / f"{up}_patches.pkl", "wb") as fh:
            pickle.dump(patches, fh)
        rows.append(ps._asdict())
    return pd.DataFrame(rows)
