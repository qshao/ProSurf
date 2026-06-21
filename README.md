# ProSurf

**Zwitterionic surface patterning as a predictor of protein thermostability**

ProSurf computes a zwitterionic surface score (Z) for protein structures and tests whether proteins with more balanced, dense, and spatially intermixed surface charges tend to be more thermostable. Validated against 7,505 human proteins with measured melting temperatures from the Meltome Atlas.

## Key result

Across 7,505 randomly sampled human proteins (AlphaFold2 v6 structures, Meltome Atlas Tm values):

| Metric | Value |
|--------|-------|
| Spearman ρ(z_mean, Tm) | **+0.165** |
| p-value | **8.1 × 10⁻⁴⁷** |
| 95% CI (bootstrap) | [0.143, 0.187] |
| Pearson r² | 2.9% |

The monotone gradient holds from 43°C to 68°C with no reversals. See [`data/prosurf_fullproteome_report.pdf`](data/prosurf_fullproteome_report.pdf) for the full analysis.

## The metric: Z = B × D̂ × M

For each surface-exposed charged residue, a local neighbourhood of radius `r = max(0.55 × Rg, 8 Å)` is evaluated:

- **B (Balance)** `= 1 − |n⁺ − n⁻| / (n⁺ + n⁻)` — penalises one-sign patches
- **D̂ (Normalised density)** `= (n⁺ + n⁻) / (π r²) / Dmax` — within-protein relative charge density
- **M (Mixing)** — distance-weighted fraction of opposite-sign neighbours within r/2

The protein-level score **z_mean** is the size-weighted mean of per-patch mean-Z after clustering high-Z locations by connected components (adjacency radius 8 Å).

## Repository structure

```
prosurf/               Python package
  config.py            MetricConfig dataclass
  io/                  AlphaFold2 fetching (EBI API)
  surface/             SASA computation, charge extraction
  metric/              Engine A (mixing_euclidean) and Engine B (cross_colocation)
  patches/             Clustering and aggregation
  validate/            AUROC, robustness sweeps, synthetic controls
  pipeline.py          analyze_structure() — single-protein entry point

scripts/
  run_all_9160.py      Full Meltome Atlas pipeline (resumable, tmux-safe)
  run_validation.py    AUROC + robustness validation on the control set

configs/
  control_set.yaml     200-protein biased control set (100 high-Tm / 100 low-Tm)

data/
  all_9160_scores.csv         Full-proteome scores (7,505 proteins)
  human_meltome_tms.csv       Meltome Atlas median Tm values
  prosurf_fullproteome_report.pdf   Full analysis report with figures
  random_1000_scores.csv      1,000-protein unbiased pilot
  pilot_200_scores.csv        200-protein biased control set scores

tests/                 pytest suite
```

## Installation

```bash
pip install -e ".[dev]"
```

Requires Python ≥ 3.11. Key dependencies: `biotite`, `numpy`, `scipy`, `requests`, `pyyaml`.

## Usage

### Score a single protein

```python
from prosurf.config import MetricConfig, PathsConfig
from prosurf.io.fetch import fetch_af2
from prosurf.pipeline import analyze_structure

cfg   = MetricConfig()
paths = PathsConfig()

pdb_path = fetch_af2("P00533", paths.data_dir)   # downloads EGFR AF2 structure
_, score = analyze_structure(pdb_path, "P00533", cfg)

print(score.z_mean, score.n_patches)
```

### Run the full Meltome Atlas pipeline

```bash
# Launches in a persistent tmux session, resumes automatically if interrupted
tmux new-session -d -s prosurf \
  'python3 scripts/run_all_9160.py 2>&1 | tee logs/run_all_9160.log'

# Monitor progress
tail -f logs/run_all_9160.log
```

### Run the biased control-set validation

```bash
python3 scripts/run_validation.py
```

### Run tests

```bash
pytest tests/
```

## Data sources

- **Thermal stability:** Meltome Atlas, Jarzab et al. 2020, *Nature Methods* 17, 495–503 — Supplementary Table S2 (55 human TPP datasets, median Tm per protein)
- **Structures:** AlphaFold2 v6, EMBL-EBI ([alphafold.ebi.ac.uk](https://alphafold.ebi.ac.uk))

## Configuration

```python
MetricConfig(
    rsasa_threshold  = 0.20,   # minimum relative SASA to be "surface"
    patch_radius_frac= 0.55,   # adaptive radius = max(frac × Rg, 8 Å)
    his_weight       = 0.0,    # histidine charge weight (0 = excluded)
    z_percentile     = 90.0,   # top-Z threshold for patch clustering
    adjacency_radius = 8.0,    # connected-component adjacency in Å
)
```
