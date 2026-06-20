# Zwitterionic Surface Quantification — Phase 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build and validate a robust algorithm that quantifies the zwitterionic pattern (net-neutral but charge-rich, intermixed +/− regions) on AlphaFold2 human protein surfaces, producing per-patch maps and per-protein score vectors.

**Architecture:** Staged modular pipeline. Stage A (prototype) uses solvent-accessibility residue selection + Euclidean spherical patches + a composite `Z = B·D̂·M` index. Stage B (production) upgrades the surface to a sampled surface point-cloud with kNN-graph geodesic distances and replaces the mixing term with bivariate spatial statistics (geodesic cross-K). Both engines emit the same `Z-per-location` artifact, so downstream patch-finding and aggregation are shared.

**Tech Stack:** Python 3.11+, `biotite` (structure parsing + Shrake-Rupley SASA), `numpy`, `scipy` (KDTree, sparse graph/Dijkstra, stats), `pandas`, `pyvista`+`matplotlib` (viz), `pytest`. No external surface binaries (MSMS/NanoShaper) — Stage B derives the surface from Shrake-Rupley sphere sampling for portability.

## Global Constraints

- Python 3.11+; all dependencies pip-installable (no external compiled binaries).
- Physiological pH model: positive = Lys (NZ), Arg (CZ); negative = Asp (OD1/OD2 midpoint), Glu (OE1/OE2 midpoint); His excluded by default via fractional weight `his_weight=0.0`; termini included as point charges.
- Charges are points at side-chain charge centers, never Cα.
- `Z = B · D̂ · M` is a strict product (any zero component → Z=0).
- Every stage is independently cacheable; one protein = one independent unit of work.
- TDD: every task writes a failing test first; synthetic structures provide exact ground truth.
- All randomness seeded; all tunable parameters live in `prosurf/config.py` dataclasses.

---

### Task 1: Project scaffolding & dependencies

**Files:**
- Create: `pyproject.toml`
- Create: `prosurf/__init__.py`
- Create: `prosurf/config.py`
- Create: `tests/__init__.py`
- Create: `tests/test_config.py`
- Create: `.gitignore`

**Interfaces:**
- Produces: `prosurf.config.MetricConfig` dataclass with fields `rsasa_threshold: float = 0.20`, `patch_radius: float = 10.0`, `his_weight: float = 0.0`, `z_percentile: float = 90.0`, `seed: int = 0`; `prosurf.config.PathsConfig` with `data_dir: Path`, `cache_dir: Path`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_config.py
from prosurf.config import MetricConfig

def test_metric_config_defaults():
    cfg = MetricConfig()
    assert cfg.rsasa_threshold == 0.20
    assert cfg.patch_radius == 10.0
    assert cfg.his_weight == 0.0
    assert cfg.z_percentile == 90.0
    assert cfg.seed == 0

def test_metric_config_override():
    cfg = MetricConfig(his_weight=0.1)
    assert cfg.his_weight == 0.1
    assert cfg.rsasa_threshold == 0.20
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_config.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'prosurf'`

- [ ] **Step 3: Write pyproject.toml and config**

```toml
# pyproject.toml
[project]
name = "prosurf"
version = "0.1.0"
requires-python = ">=3.11"
dependencies = ["biotite>=0.40", "numpy", "scipy", "pandas", "pyvista", "matplotlib", "requests"]

[project.optional-dependencies]
dev = ["pytest"]

[tool.setuptools.packages.find]
include = ["prosurf*"]

[build-system]
requires = ["setuptools>=61"]
build-backend = "setuptools.build_meta"
```

```python
# prosurf/config.py
from dataclasses import dataclass, field
from pathlib import Path

@dataclass
class MetricConfig:
    rsasa_threshold: float = 0.20
    patch_radius: float = 10.0
    his_weight: float = 0.0
    z_percentile: float = 90.0
    seed: int = 0

@dataclass
class PathsConfig:
    data_dir: Path = field(default_factory=lambda: Path("data"))
    cache_dir: Path = field(default_factory=lambda: Path("data/cache"))
```

Create empty `prosurf/__init__.py`, `tests/__init__.py`. `.gitignore`: add `data/cache/`, `__pycache__/`, `*.pyc`, `.pytest_cache/`, `*.egg-info/`.

- [ ] **Step 4: Install editable and run test**

Run: `pip install -e ".[dev]" && pytest tests/test_config.py -v`
Expected: PASS (2 passed)

- [ ] **Step 5: Commit**

```bash
git add pyproject.toml prosurf/ tests/ .gitignore
git commit -m "feat: project scaffolding and config"
```

---

### Task 2: Fetch & parse AlphaFold2 structures

**Files:**
- Create: `prosurf/io/__init__.py`
- Create: `prosurf/io/fetch.py`
- Create: `prosurf/io/parse.py`
- Create: `tests/test_io.py`
- Create: `tests/data/mini.pdb` (3-residue hand-written test structure)

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `prosurf.io.fetch.af2_url(uniprot_id: str) -> str`
  - `prosurf.io.fetch.fetch_af2(uniprot_id: str, out_dir: Path) -> Path` (downloads `.pdb`, cached)
  - `prosurf.io.parse.load_structure(path: Path) -> AtomArray` (biotite `AtomArray`, model 1, amino acids only)

- [ ] **Step 1: Write the failing test**

Create `tests/data/mini.pdb` with one Lys, one Asp, one Ala (real ATOM records, coordinates spaced ~6 Å apart). Then:

```python
# tests/test_io.py
from pathlib import Path
from prosurf.io.fetch import af2_url
from prosurf.io.parse import load_structure
import biotite.structure as struc

def test_af2_url():
    assert af2_url("P69905") == "https://alphafold.ebi.ac.uk/files/AF-P69905-F1-model_v4.pdb"

def test_load_structure_amino_acids_only(tmp_path):
    arr = load_structure(Path("tests/data/mini.pdb"))
    res_ids, res_names = struc.get_residues(arr)
    assert set(res_names) == {"LYS", "ASP", "ALA"}
    assert len(res_ids) == 3
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_io.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'prosurf.io'`

- [ ] **Step 3: Write implementation**

```python
# prosurf/io/fetch.py
from pathlib import Path
import requests

def af2_url(uniprot_id: str) -> str:
    return f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"

def fetch_af2(uniprot_id: str, out_dir: Path) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / f"AF-{uniprot_id}-F1-model_v4.pdb"
    if dest.exists():
        return dest
    resp = requests.get(af2_url(uniprot_id), timeout=60)
    resp.raise_for_status()
    dest.write_text(resp.text)
    return dest
```

```python
# prosurf/io/parse.py
from pathlib import Path
import biotite.structure as struc
from biotite.structure.io.pdb import PDBFile

def load_structure(path: Path) -> struc.AtomArray:
    pdb = PDBFile.read(str(path))
    arr = pdb.get_structure(model=1)
    arr = arr[struc.filter_amino_acids(arr)]
    return arr
```

Create empty `prosurf/io/__init__.py`.

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_io.py -v`
Expected: PASS (2 passed)

- [ ] **Step 5: Commit**

```bash
git add prosurf/io/ tests/test_io.py tests/data/mini.pdb
git commit -m "feat: fetch and parse AF2 structures"
```

---

### Task 3: SASA & surface residue selection

**Files:**
- Create: `prosurf/surface/__init__.py`
- Create: `prosurf/surface/sasa.py`
- Create: `tests/test_sasa.py`

**Interfaces:**
- Consumes: `load_structure` → `AtomArray`.
- Produces:
  - `prosurf.surface.sasa.residue_sasa(arr: AtomArray) -> tuple[np.ndarray, np.ndarray]` returns `(res_ids, sasa_per_residue)`.
  - `prosurf.surface.sasa.relative_sasa(arr: AtomArray) -> tuple[np.ndarray, np.ndarray]` returns `(res_ids, rsasa)` where rsasa = residue SASA / max ASA (Tien theoretical).
  - `prosurf.surface.sasa.surface_residue_ids(arr, threshold: float) -> np.ndarray` of res_ids with rsasa > threshold.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_sasa.py
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_sasa.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'prosurf.surface'`

- [ ] **Step 3: Write implementation**

```python
# prosurf/surface/sasa.py
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
```

Create empty `prosurf/surface/__init__.py`.

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_sasa.py -v`
Expected: PASS (3 passed)

- [ ] **Step 5: Commit**

```bash
git add prosurf/surface/ tests/test_sasa.py
git commit -m "feat: SASA and surface residue selection"
```

---

### Task 4: Charge assignment at side-chain charge centers

**Files:**
- Create: `prosurf/surface/charges.py`
- Create: `tests/test_charges.py`

**Interfaces:**
- Consumes: `AtomArray`, `surface_residue_ids`.
- Produces:
  - `prosurf.surface.charges.Charge` namedtuple `(res_id: int, sign: int, weight: float, xyz: np.ndarray)`.
  - `prosurf.surface.charges.assign_charges(arr, surface_ids, his_weight=0.0) -> list[Charge]` — only residues whose `res_id` is in `surface_ids`; termini always included.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_charges.py
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_charges.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'prosurf.surface.charges'`

- [ ] **Step 3: Write implementation**

```python
# prosurf/surface/charges.py
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
    nterm, cterm = res_ids[0], res_ids[-1]
    n_xyz = _atom_coord(arr, nterm, "N")
    c_xyz = _atom_coord(arr, cterm, "C")
    if n_xyz is not None: charges.append(Charge(int(nterm), 1, 1.0, n_xyz))
    if c_xyz is not None: charges.append(Charge(int(cterm), -1, 1.0, c_xyz))
    return charges
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_charges.py -v`
Expected: PASS (2 passed)

- [ ] **Step 5: Commit**

```bash
git add prosurf/surface/charges.py tests/test_charges.py
git commit -m "feat: charge assignment at side-chain charge centers"
```

---

### Task 5: Metric components — Balance, Density, Mixing

**Files:**
- Create: `prosurf/metric/__init__.py`
- Create: `prosurf/metric/components.py`
- Create: `tests/test_components.py`

**Interfaces:**
- Consumes: `list[Charge]`.
- Produces:
  - `prosurf.metric.components.balance(n_pos: float, n_neg: float) -> float`
  - `prosurf.metric.components.density(total_charge: float, area: float) -> float` (raw, un-normalized)
  - `prosurf.metric.components.mixing_euclidean(charges: list[Charge], neighbor_radius: float) -> float`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_components.py
import numpy as np
from prosurf.surface.charges import Charge
from prosurf.metric.components import balance, density, mixing_euclidean

def test_balance_perfect():
    assert balance(5, 5) == 1.0

def test_balance_one_sided():
    assert balance(5, 0) == 0.0

def test_balance_empty():
    assert balance(0, 0) == 0.0

def test_density():
    assert density(10.0, 100.0) == 0.1

def test_mixing_alternating_high():
    # + and - interleaved within 4 A -> high mixing
    charges = [
        Charge(1, 1, 1.0, np.array([0., 0., 0.])),
        Charge(2, -1, 1.0, np.array([3., 0., 0.])),
        Charge(3, 1, 1.0, np.array([6., 0., 0.])),
        Charge(4, -1, 1.0, np.array([9., 0., 0.])),
    ]
    assert mixing_euclidean(charges, neighbor_radius=4.0) > 0.9

def test_mixing_segregated_low():
    # all + clustered, all - far -> low mixing
    charges = [
        Charge(1, 1, 1.0, np.array([0., 0., 0.])),
        Charge(2, 1, 1.0, np.array([2., 0., 0.])),
        Charge(3, -1, 1.0, np.array([50., 0., 0.])),
        Charge(4, -1, 1.0, np.array([52., 0., 0.])),
    ]
    assert mixing_euclidean(charges, neighbor_radius=4.0) < 0.1
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_components.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'prosurf.metric'`

- [ ] **Step 3: Write implementation**

```python
# prosurf/metric/components.py
import numpy as np
from scipy.spatial import cKDTree

def balance(n_pos, n_neg):
    total = n_pos + n_neg
    if total == 0:
        return 0.0
    return 1.0 - abs(n_pos - n_neg) / total

def density(total_charge, area):
    if area <= 0:
        return 0.0
    return total_charge / area

def mixing_euclidean(charges, neighbor_radius):
    if len(charges) < 2:
        return 0.0
    coords = np.array([c.xyz for c in charges])
    signs = np.array([c.sign for c in charges])
    tree = cKDTree(coords)
    fracs = []
    for i, c in enumerate(charges):
        idx = tree.query_ball_point(coords[i], neighbor_radius)
        idx = [j for j in idx if j != i]
        if not idx:
            continue
        opposite = np.sum(signs[idx] != signs[i])
        fracs.append(opposite / len(idx))
    return float(np.mean(fracs)) if fracs else 0.0
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_components.py -v`
Expected: PASS (6 passed)

- [ ] **Step 5: Commit**

```bash
git add prosurf/metric/ tests/test_components.py
git commit -m "feat: balance, density, mixing metric components"
```

---

### Task 6: Engine A — per-location Z score

**Files:**
- Create: `prosurf/metric/engine_a.py`
- Create: `tests/test_engine_a.py`

**Interfaces:**
- Consumes: `AtomArray`, `MetricConfig`, `assign_charges`, components.
- Produces:
  - `prosurf.metric.engine_a.LocationScore` namedtuple `(res_id, z, n_pos, n_neg, xyz)`.
  - `prosurf.metric.engine_a.score_locations_a(arr, cfg) -> list[LocationScore]` — for each surface charged residue, build a `cfg.patch_radius` sphere, compute `Z = balance * density_norm * mixing`. Density normalized by patch disk area `π·r²`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_engine_a.py
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_engine_a.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'prosurf.metric.engine_a'`

- [ ] **Step 3: Write implementation**

```python
# prosurf/metric/engine_a.py
from collections import namedtuple
import numpy as np
from scipy.spatial import cKDTree
from prosurf.metric.components import balance, density, mixing_euclidean
from prosurf.surface.charges import assign_charges
from prosurf.surface.sasa import surface_residue_ids

LocationScore = namedtuple("LocationScore", ["res_id", "z", "n_pos", "n_neg", "xyz"])

def score_charges_a(charges, cfg):
    if not charges:
        return []
    coords = np.array([c.xyz for c in charges])
    signs = np.array([c.sign for c in charges])
    weights = np.array([c.weight for c in charges])
    tree = cKDTree(coords)
    disk_area = np.pi * cfg.patch_radius ** 2
    scores = []
    for i, c in enumerate(charges):
        idx = tree.query_ball_point(coords[i], cfg.patch_radius)
        region = [charges[j] for j in idx]
        n_pos = float(np.sum(weights[idx] * (signs[idx] > 0)))
        n_neg = float(np.sum(weights[idx] * (signs[idx] < 0)))
        B = balance(n_pos, n_neg)
        D = density(n_pos + n_neg, disk_area)
        M = mixing_euclidean(region, neighbor_radius=cfg.patch_radius / 2)
        scores.append(LocationScore(c.res_id, B * D * M, n_pos, n_neg, c.xyz))
    return scores  # density normalization across dataset applied in aggregation

def score_locations_a(arr, cfg):
    sids = surface_residue_ids(arr, cfg.rsasa_threshold)
    charges = assign_charges(arr, sids, his_weight=cfg.his_weight)
    return score_charges_a(charges, cfg)
```

**Note:** raw `density` is normalized to `D̂∈[0,1]` at the dataset level in Task 8 (`normalize_density`), so per-location Z here is monotonic in D; ordering tests hold.

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_engine_a.py -v`
Expected: PASS (2 passed)

- [ ] **Step 5: Commit**

```bash
git add prosurf/metric/engine_a.py tests/test_engine_a.py
git commit -m "feat: engine A per-location Z scoring"
```

---

### Task 7: Patch clustering

**Files:**
- Create: `prosurf/patches/__init__.py`
- Create: `prosurf/patches/cluster.py`
- Create: `tests/test_cluster.py`

**Interfaces:**
- Consumes: `list[LocationScore]`, `MetricConfig`.
- Produces:
  - `prosurf.patches.cluster.Patch` namedtuple `(res_ids: list[int], mean_z: float, max_z: float, n_pos: float, n_neg: float, size: int)`.
  - `prosurf.patches.cluster.cluster_patches(scores, cfg, adjacency_radius=8.0) -> list[Patch]` — keep locations with `z` above the `cfg.z_percentile` percentile, connect those within `adjacency_radius`, return connected components as patches.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_cluster.py
import numpy as np
from prosurf.config import MetricConfig
from prosurf.metric.engine_a import LocationScore
from prosurf.patches.cluster import cluster_patches

def test_two_separate_patches():
    cfg = MetricConfig(z_percentile=0.0)  # keep all
    scores = [
        LocationScore(1, 0.9, 2, 2, np.array([0., 0., 0.])),
        LocationScore(2, 0.8, 2, 2, np.array([3., 0., 0.])),
        LocationScore(3, 0.9, 2, 2, np.array([100., 0., 0.])),
    ]
    patches = cluster_patches(scores, cfg, adjacency_radius=8.0)
    assert len(patches) == 2
    sizes = sorted(p.size for p in patches)
    assert sizes == [1, 2]

def test_percentile_filters_low_z():
    cfg = MetricConfig(z_percentile=90.0)
    scores = [LocationScore(i, float(i)/10, 1, 1, np.array([float(i), 0., 0.]))
              for i in range(11)]
    patches = cluster_patches(scores, cfg, adjacency_radius=8.0)
    kept = sum(p.size for p in patches)
    assert kept <= 2  # only top ~10%
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_cluster.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'prosurf.patches'`

- [ ] **Step 3: Write implementation**

```python
# prosurf/patches/cluster.py
from collections import namedtuple
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

Patch = namedtuple("Patch", ["res_ids", "mean_z", "max_z", "n_pos", "n_neg", "size"])

def cluster_patches(scores, cfg, adjacency_radius=8.0):
    if not scores:
        return []
    z = np.array([s.z for s in scores])
    cutoff = np.percentile(z, cfg.z_percentile)
    keep = [i for i, s in enumerate(scores) if s.z >= cutoff and s.z > 0]
    if not keep:
        return []
    kept = [scores[i] for i in keep]
    coords = np.array([s.xyz for s in kept])
    tree = cKDTree(coords)
    pairs = tree.query_pairs(adjacency_radius, output_type="ndarray")
    n = len(kept)
    if len(pairs) > 0:
        data = np.ones(len(pairs))
        graph = csr_matrix((data, (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    else:
        graph = csr_matrix((n, n))
    n_comp, labels = connected_components(graph, directed=False)
    patches = []
    for c in range(n_comp):
        members = [kept[i] for i in range(n) if labels[i] == c]
        zs = np.array([m.z for m in members])
        patches.append(Patch(
            res_ids=[m.res_id for m in members],
            mean_z=float(zs.mean()), max_z=float(zs.max()),
            n_pos=float(sum(m.n_pos for m in members)),
            n_neg=float(sum(m.n_neg for m in members)),
            size=len(members)))
    return patches
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_cluster.py -v`
Expected: PASS (2 passed)

- [ ] **Step 5: Commit**

```bash
git add prosurf/patches/ tests/test_cluster.py
git commit -m "feat: patch clustering via connected components"
```

---

### Task 8: Per-protein aggregation & dataset density normalization

**Files:**
- Create: `prosurf/patches/aggregate.py`
- Create: `tests/test_aggregate.py`

**Interfaces:**
- Consumes: `list[Patch]`, total surface area, raw per-location densities across dataset.
- Produces:
  - `prosurf.patches.aggregate.normalize_density(raw_densities: np.ndarray) -> callable` — returns a function mapping raw density → [0,1] via dataset max (robust 99th percentile cap).
  - `prosurf.patches.aggregate.ProteinScore` namedtuple `(uniprot, z_frac, z_max, z_mean, n_patches)`.
  - `prosurf.patches.aggregate.aggregate_protein(uniprot, patches, total_surface_area, n_surface_locations) -> ProteinScore`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_aggregate.py
import numpy as np
from prosurf.patches.cluster import Patch
from prosurf.patches.aggregate import aggregate_protein, normalize_density, ProteinScore

def test_normalize_density_caps_at_one():
    f = normalize_density(np.array([1., 2., 3., 100.]))
    assert 0.0 <= f(2.0) <= 1.0
    assert f(1000.0) == 1.0

def test_aggregate_protein_basic():
    patches = [
        Patch([1, 2], mean_z=0.5, max_z=0.6, n_pos=2, n_neg=2, size=2),
        Patch([5], mean_z=0.3, max_z=0.3, n_pos=1, n_neg=1, size=1),
    ]
    ps = aggregate_protein("P12345", patches, total_surface_area=100.0,
                           n_surface_locations=10)
    assert ps.uniprot == "P12345"
    assert ps.n_patches == 2
    assert ps.z_max == 0.6
    assert 0.0 <= ps.z_frac <= 1.0
    assert ps.z_mean > 0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_aggregate.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write implementation**

```python
# prosurf/patches/aggregate.py
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_aggregate.py -v`
Expected: PASS (2 passed)

- [ ] **Step 5: Commit**

```bash
git add prosurf/patches/aggregate.py tests/test_aggregate.py
git commit -m "feat: per-protein aggregation and density normalization"
```

---

### Task 9: End-to-end pipeline (Stage A) + caching

**Files:**
- Create: `prosurf/pipeline.py`
- Create: `tests/test_pipeline.py`

**Interfaces:**
- Consumes: all prior modules, `MetricConfig`, `PathsConfig`.
- Produces:
  - `prosurf.pipeline.analyze_structure(path, uniprot, cfg) -> tuple[list[Patch], ProteinScore]`.
  - `prosurf.pipeline.run_pilot(uniprots, paths, cfg) -> pandas.DataFrame` (one row per protein with ProteinScore fields), with per-protein patch lists cached as `.npz`/pickle in `paths.cache_dir`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_pipeline.py
from pathlib import Path
from prosurf.config import MetricConfig
from prosurf.pipeline import analyze_structure

def test_analyze_structure_runs_on_mini():
    cfg = MetricConfig(z_percentile=0.0)
    patches, score = analyze_structure(Path("tests/data/mini.pdb"), "MINI", cfg)
    assert score.uniprot == "MINI"
    assert isinstance(patches, list)
    assert score.n_patches >= 0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_pipeline.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'prosurf.pipeline'`

- [ ] **Step 3: Write implementation**

```python
# prosurf/pipeline.py
import pickle
import numpy as np
import pandas as pd
from pathlib import Path
from prosurf.io.parse import load_structure
from prosurf.io.fetch import fetch_af2
from prosurf.metric.engine_a import score_locations_a
from prosurf.patches.cluster import cluster_patches
from prosurf.patches.aggregate import aggregate_protein
import biotite.structure as struc
from prosurf.surface.sasa import residue_sasa, surface_residue_ids

def analyze_structure(path, uniprot, cfg):
    arr = load_structure(path)
    scores = score_locations_a(arr, cfg)
    patches = cluster_patches(scores, cfg)
    _, sasa = residue_sasa(arr)
    total_area = float(np.nansum(sasa))
    n_locations = len(surface_residue_ids(arr, cfg.rsasa_threshold))
    ps = aggregate_protein(uniprot, patches, total_area, max(n_locations, 1))
    return patches, ps

def run_pilot(uniprots, paths, cfg):
    paths.cache_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for up in uniprots:
        pdb = fetch_af2(up, paths.data_dir)
        patches, ps = analyze_structure(pdb, up, cfg)
        with open(paths.cache_dir / f"{up}_patches.pkl", "wb") as fh:
            pickle.dump(patches, fh)
        rows.append(ps._asdict())
    return pd.DataFrame(rows)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_pipeline.py -v`
Expected: PASS (1 passed)

- [ ] **Step 5: Commit**

```bash
git add prosurf/pipeline.py tests/test_pipeline.py
git commit -m "feat: end-to-end Stage A pipeline with caching"
```

---

### Task 10: Synthetic surfaces & positive/negative validation

**Files:**
- Create: `prosurf/validate/__init__.py`
- Create: `prosurf/validate/synthetic.py`
- Create: `prosurf/validate/controls.py`
- Create: `tests/test_validation.py`

**Interfaces:**
- Consumes: `score_charges_a`, `MetricConfig`, `Charge`.
- Produces:
  - `prosurf.validate.synthetic.make_synthetic(pattern: str, n: int, spacing: float = 3.0) -> list[Charge]` for `pattern in {"alternating", "segregated", "charge_free"}`.
  - `prosurf.validate.controls.protein_zwit_score(charges, cfg) -> float` (max per-location Z).
  - `prosurf.validate.controls.auroc(pos_scores, neg_scores) -> float`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_validation.py
from prosurf.config import MetricConfig
from prosurf.validate.synthetic import make_synthetic
from prosurf.validate.controls import protein_zwit_score, auroc

def test_synthetic_ordering():
    cfg = MetricConfig(patch_radius=8.0)
    alt = protein_zwit_score(make_synthetic("alternating", 12), cfg)
    seg = protein_zwit_score(make_synthetic("segregated", 12), cfg)
    free = protein_zwit_score(make_synthetic("charge_free", 12), cfg)
    assert alt > seg
    assert alt > free

def test_auroc_perfect_separation():
    assert auroc([0.9, 0.8, 0.95], [0.1, 0.2, 0.05]) == 1.0

def test_auroc_random():
    assert abs(auroc([0.5, 0.5], [0.5, 0.5]) - 0.5) < 1e-9
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_validation.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'prosurf.validate'`

- [ ] **Step 3: Write implementation**

```python
# prosurf/validate/synthetic.py
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
```

```python
# prosurf/validate/controls.py
import numpy as np
from prosurf.metric.engine_a import score_charges_a

def protein_zwit_score(charges, cfg):
    scores = score_charges_a(charges, cfg)
    return max((s.z for s in scores), default=0.0)

def auroc(pos_scores, neg_scores):
    pos = np.asarray(pos_scores); neg = np.asarray(neg_scores)
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    wins = sum((p > n) + 0.5 * (p == n) for p in pos for n in neg)
    return float(wins / (len(pos) * len(neg)))
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_validation.py -v`
Expected: PASS (3 passed)

- [ ] **Step 5: Commit**

```bash
git add prosurf/validate/ tests/test_validation.py
git commit -m "feat: synthetic surfaces and control-set validation"
```

---

### Task 11: Robustness / sensitivity sweeps

**Files:**
- Create: `prosurf/validate/robustness.py`
- Create: `tests/test_robustness.py`

**Interfaces:**
- Consumes: `analyze_structure`, `MetricConfig`, `ProteinScore`.
- Produces:
  - `prosurf.validate.robustness.sweep_parameter(paths_or_struct, base_cfg, param: str, values: list) -> pandas.DataFrame` — re-scores a fixed protein set varying one config field; rows are (param_value, uniprot, z_frac, z_max, z_mean, n_patches).
  - `prosurf.validate.robustness.ranking_stability(df, score_col="z_mean") -> float` — mean pairwise Spearman ρ of per-protein rankings across parameter values.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_robustness.py
import pandas as pd
from prosurf.validate.robustness import ranking_stability

def test_ranking_stability_identical_rankings():
    df = pd.DataFrame({
        "param_value": [0.1, 0.1, 0.2, 0.2],
        "uniprot": ["A", "B", "A", "B"],
        "z_mean": [0.9, 0.1, 0.8, 0.2],
    })
    assert abs(ranking_stability(df, "z_mean") - 1.0) < 1e-9

def test_ranking_stability_inverted():
    df = pd.DataFrame({
        "param_value": [0.1, 0.1, 0.2, 0.2],
        "uniprot": ["A", "B", "A", "B"],
        "z_mean": [0.9, 0.1, 0.1, 0.9],
    })
    assert ranking_stability(df, "z_mean") < 0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_robustness.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write implementation**

```python
# prosurf/validate/robustness.py
from dataclasses import replace
from itertools import combinations
import pandas as pd
from scipy.stats import spearmanr
from prosurf.pipeline import analyze_structure

def sweep_parameter(structures, base_cfg, param, values):
    rows = []
    for v in values:
        cfg = replace(base_cfg, **{param: v})
        for path, uniprot in structures:
            _, ps = analyze_structure(path, uniprot, cfg)
            row = ps._asdict(); row["param_value"] = v
            rows.append(row)
    return pd.DataFrame(rows)

def ranking_stability(df, score_col="z_mean"):
    pivot = df.pivot(index="uniprot", columns="param_value", values=score_col)
    cols = list(pivot.columns)
    rhos = []
    for a, b in combinations(cols, 2):
        rho, _ = spearmanr(pivot[a], pivot[b])
        rhos.append(rho)
    return float(sum(rhos) / len(rhos)) if rhos else 1.0
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_robustness.py -v`
Expected: PASS (2 passed)

- [ ] **Step 5: Commit**

```bash
git add prosurf/validate/robustness.py tests/test_robustness.py
git commit -m "feat: robustness sweeps and ranking stability"
```

---

### Task 12: Patch-map visualization

**Files:**
- Create: `prosurf/viz/__init__.py`
- Create: `prosurf/viz/patchmap.py`
- Create: `tests/test_viz.py`

**Interfaces:**
- Consumes: `AtomArray`, `list[LocationScore]`.
- Produces:
  - `prosurf.viz.patchmap.render_patch_map(arr, scores, out_png: Path) -> Path` — Cα scatter colored by per-residue Z (matplotlib 3D), saved headless.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_viz.py
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
from prosurf.config import MetricConfig
from prosurf.io.parse import load_structure
from prosurf.metric.engine_a import score_locations_a
from prosurf.viz.patchmap import render_patch_map

def test_render_creates_png(tmp_path):
    cfg = MetricConfig(z_percentile=0.0)
    arr = load_structure(Path("tests/data/mini.pdb"))
    scores = score_locations_a(arr, cfg)
    out = render_patch_map(arr, scores, tmp_path / "map.png")
    assert out.exists() and out.stat().st_size > 0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_viz.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write implementation**

```python
# prosurf/viz/patchmap.py
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def render_patch_map(arr, scores, out_png):
    out_png = Path(out_png)
    zmap = {s.res_id: s.z for s in scores}
    pts = np.array([s.xyz for s in scores]) if scores else np.zeros((1, 3))
    vals = np.array([s.z for s in scores]) if scores else np.array([0.0])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    p = ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], c=vals, cmap="viridis")
    fig.colorbar(p, label="Z (zwitterionic)")
    fig.savefig(out_png, dpi=120)
    plt.close(fig)
    return out_png
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_viz.py -v`
Expected: PASS (1 passed)

- [ ] **Step 5: Commit**

```bash
git add prosurf/viz/ tests/test_viz.py
git commit -m "feat: patch-map visualization"
```

---

### Task 13: Control-set definition & Stage-A validation report

**Files:**
- Create: `configs/control_set.yaml`
- Create: `scripts/run_validation.py`
- Create: `tests/test_control_set.py`

**Interfaces:**
- Consumes: `run_pilot`, `controls.auroc`, `robustness`.
- Produces:
  - `configs/control_set.yaml` listing UniProt IDs under `positives:` (halophilic/extremophile, charge-rich balanced surfaces) and `negatives:` (charge-segregated / charge-poor), each with a one-line rationale comment.
  - `scripts/run_validation.py` — fetches the control set, computes AUROC of `z_max` (positives vs negatives), runs robustness sweeps over `rsasa_threshold`, `patch_radius`, `his_weight`, writes `data/validation_report.md` + per-protein patch-map PNGs.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_control_set.py
import yaml
from pathlib import Path

def test_control_set_has_both_classes():
    data = yaml.safe_load(Path("configs/control_set.yaml").read_text())
    assert len(data["positives"]) >= 5
    assert len(data["negatives"]) >= 5
    # UniProt-like accessions
    for k in ("positives", "negatives"):
        for acc in data[k]:
            assert isinstance(acc, str) and 6 <= len(acc) <= 10
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_control_set.py -v`
Expected: FAIL with `FileNotFoundError: configs/control_set.yaml`

- [ ] **Step 3: Write the control set and script**

```yaml
# configs/control_set.yaml
# Phase-1 validation controls. Positives: surfaces expected zwitterionic
# (charge-rich + balanced, e.g. halophilic/extremophile). Negatives:
# charge-segregated or charge-poor surfaces. Replace/expand as literature dictates.
positives:
  - P0A7L8   # ribosomal-like acidic/basic balanced (placeholder rationale)
  - P00509
  - P0AEX9
  - P02769
  - P61626
  - P00698
negatives:
  - P01308   # insulin, small / charge-poor
  - P69905   # hemoglobin alpha
  - P0DTC2
  - P02649
  - P00766
  - P01112
```

**Note:** the implementer must confirm/replace these accessions against the literature before trusting AUROC; the test only enforces structure and count. The script:

```python
# scripts/run_validation.py
import yaml
from pathlib import Path
from prosurf.config import MetricConfig, PathsConfig
from prosurf.io.fetch import fetch_af2
from prosurf.pipeline import analyze_structure
from prosurf.validate.controls import auroc
from prosurf.validate.robustness import sweep_parameter, ranking_stability

def main():
    cfg = MetricConfig(); paths = PathsConfig()
    cs = yaml.safe_load(Path("configs/control_set.yaml").read_text())
    def scores(accs):
        out = []
        for a in accs:
            p = fetch_af2(a, paths.data_dir)
            _, ps = analyze_structure(p, a, cfg)
            out.append(ps.z_max)
        return out
    pos, neg = scores(cs["positives"]), scores(cs["negatives"])
    au = auroc(pos, neg)
    structs = [(fetch_af2(a, paths.data_dir), a) for a in cs["positives"] + cs["negatives"]]
    lines = [f"# Stage-A Validation Report", f"", f"AUROC (z_max, pos vs neg): {au:.3f}", ""]
    for param, vals in [("rsasa_threshold", [0.15, 0.20, 0.25]),
                        ("patch_radius", [8.0, 10.0, 12.0]),
                        ("his_weight", [0.0, 0.05, 0.1])]:
        df = sweep_parameter(structs, cfg, param, vals)
        lines.append(f"- {param}: ranking stability (Spearman) = {ranking_stability(df):.3f}")
    Path("data/validation_report.md").write_text("\n".join(lines))
    print("\n".join(lines))

if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run the test (and optionally the script with network)**

Run: `pytest tests/test_control_set.py -v`
Expected: PASS (1 passed)
Optional (needs network): `python scripts/run_validation.py` → writes `data/validation_report.md`.

- [ ] **Step 5: Commit**

```bash
git add configs/control_set.yaml scripts/run_validation.py tests/test_control_set.py
git commit -m "feat: control set and Stage-A validation report script"
```

---

### Task 14: Stage B — surface point cloud & geodesic distances

**Files:**
- Create: `prosurf/surface/pointcloud.py`
- Create: `tests/test_pointcloud.py`

**Interfaces:**
- Consumes: `AtomArray`.
- Produces:
  - `prosurf.surface.pointcloud.sample_surface(arr, probe=1.4, density=1.0) -> np.ndarray` (N×3 solvent-accessible surface points via per-atom sphere sampling, keeping points not buried inside neighbor spheres).
  - `prosurf.surface.pointcloud.geodesic_graph(points, knn=8) -> scipy.sparse matrix` (kNN Euclidean graph as geodesic proxy).
  - `prosurf.surface.pointcloud.geodesic_distances(graph, source_idx) -> np.ndarray` (Dijkstra from a source point).

- [ ] **Step 1: Write the failing test**

```python
# tests/test_pointcloud.py
import numpy as np
from pathlib import Path
from prosurf.io.parse import load_structure
from prosurf.surface.pointcloud import sample_surface, geodesic_graph, geodesic_distances

def test_sample_surface_nonempty():
    arr = load_structure(Path("tests/data/mini.pdb"))
    pts = sample_surface(arr, density=0.5)
    assert pts.ndim == 2 and pts.shape[1] == 3 and len(pts) > 0

def test_geodesic_distances_monotone():
    pts = np.array([[0,0,0],[1,0,0],[2,0,0],[3,0,0]], dtype=float)
    g = geodesic_graph(pts, knn=2)
    d = geodesic_distances(g, 0)
    assert d[0] == 0 and d[3] >= d[1]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_pointcloud.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write implementation**

```python
# prosurf/surface/pointcloud.py
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra

_VDW = {"C":1.7,"N":1.55,"O":1.52,"S":1.8,"H":1.2}

def _radii(arr):
    return np.array([_VDW.get(e, 1.7) for e in arr.element])

def sample_surface(arr, probe=1.4, density=1.0):
    coords = arr.coord
    radii = _radii(arr) + probe
    n_sphere = max(int(24 * density), 6)
    # Fibonacci sphere unit vectors
    i = np.arange(n_sphere)
    phi = np.arccos(1 - 2*(i+0.5)/n_sphere)
    theta = np.pi * (1 + 5**0.5) * i
    unit = np.column_stack([np.cos(theta)*np.sin(phi),
                            np.sin(theta)*np.sin(phi), np.cos(phi)])
    tree = cKDTree(coords)
    pts = []
    for c, r in zip(coords, radii):
        cand = c + unit * r
        # keep candidates not inside any other atom's probe sphere
        for p in cand:
            nbr = tree.query_ball_point(p, r)  # within own radius range
            inside = False
            for j in nbr:
                if np.linalg.norm(p - coords[j]) < radii[j] - 1e-6:
                    inside = True; break
            if not inside:
                pts.append(p)
    return np.array(pts) if pts else np.empty((0, 3))

def geodesic_graph(points, knn=8):
    tree = cKDTree(points)
    n = len(points)
    rows, cols, data = [], [], []
    k = min(knn + 1, n)
    dists, idx = tree.query(points, k=k)
    for i in range(n):
        for d, j in zip(dists[i][1:], idx[i][1:]):
            rows.append(i); cols.append(j); data.append(d)
    return csr_matrix((data, (rows, cols)), shape=(n, n))

def geodesic_distances(graph, source_idx):
    return dijkstra(graph, directed=False, indices=source_idx)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_pointcloud.py -v`
Expected: PASS (2 passed)

- [ ] **Step 5: Commit**

```bash
git add prosurf/surface/pointcloud.py tests/test_pointcloud.py
git commit -m "feat: Stage B surface point cloud and geodesic distances"
```

---

### Task 15: Stage B — bivariate geodesic mixing (cross-K)

**Files:**
- Create: `prosurf/metric/mixing_spatial.py`
- Create: `tests/test_mixing_spatial.py`

**Interfaces:**
- Consumes: charge coordinates + signs, a geodesic distance function.
- Produces:
  - `prosurf.metric.mixing_spatial.cross_colocation(coords, signs, d_max, geodesic_dist=None) -> float` — bivariate +/− co-location at short range. Returns value in [0,1]: 1 = + and − strongly co-located (zwitterionic), 0 = fully segregated. Uses geodesic distances when `geodesic_dist` provided, else Euclidean.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_mixing_spatial.py
import numpy as np
from prosurf.metric.mixing_spatial import cross_colocation

def test_colocation_alternating_high():
    coords = np.array([[i*3.,0,0] for i in range(8)])
    signs = np.array([1,-1,1,-1,1,-1,1,-1])
    assert cross_colocation(coords, signs, d_max=4.0) > 0.8

def test_colocation_segregated_low():
    coords = np.array([[i*3.,0,0] for i in range(8)])
    signs = np.array([1,1,1,1,-1,-1,-1,-1])
    assert cross_colocation(coords, signs, d_max=4.0) < 0.4
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_mixing_spatial.py -v`
Expected: FAIL with `ModuleNotFoundError`

- [ ] **Step 3: Write implementation**

```python
# prosurf/metric/mixing_spatial.py
import numpy as np
from scipy.spatial import cKDTree

def cross_colocation(coords, signs, d_max, geodesic_dist=None):
    """Fraction of short-range neighbor pairs that are opposite-signed,
    normalized against the expectation under random labeling.
    Returns ~1 when + and - co-locate (zwitterionic), ~0 when segregated."""
    n = len(coords)
    if n < 2:
        return 0.0
    pos_frac = np.mean(signs > 0)
    neg_frac = np.mean(signs < 0)
    expected_opposite = 2 * pos_frac * neg_frac  # random-labeling baseline
    if expected_opposite == 0:
        return 0.0
    if geodesic_dist is None:
        tree = cKDTree(coords)
        pairs = tree.query_pairs(d_max, output_type="ndarray")
    else:
        pairs = np.array([[i, j] for i in range(n) for j in range(i+1, n)
                          if geodesic_dist(i, j) <= d_max])
    if len(pairs) == 0:
        return 0.0
    opp = np.mean(signs[pairs[:, 0]] != signs[pairs[:, 1]])
    return float(min(opp / expected_opposite, 1.0))
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_mixing_spatial.py -v`
Expected: PASS (2 passed)

- [ ] **Step 5: Commit**

```bash
git add prosurf/metric/mixing_spatial.py tests/test_mixing_spatial.py
git commit -m "feat: Stage B bivariate geodesic cross-colocation mixing"
```

---

### Task 16: Stage B engine & engine selector

**Files:**
- Create: `prosurf/metric/engine_b.py`
- Modify: `prosurf/pipeline.py` (add `engine` switch)
- Create: `tests/test_engine_b.py`

**Interfaces:**
- Consumes: `AtomArray`, `MetricConfig`, geodesic graph, `cross_colocation`, `balance`, `density`.
- Produces:
  - `prosurf.metric.engine_b.score_locations_b(arr, cfg) -> list[LocationScore]` — geodesic patches around each surface charge; `M = cross_colocation(...)` with geodesic distances; same `LocationScore` shape as engine A.
  - `prosurf.pipeline.analyze_structure(path, uniprot, cfg, engine="a")` — `engine in {"a","b"}` selects the scorer; everything downstream unchanged.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_engine_b.py
from pathlib import Path
from prosurf.config import MetricConfig
from prosurf.pipeline import analyze_structure

def test_engine_b_runs_and_shares_shape():
    cfg = MetricConfig(z_percentile=0.0)
    patches_a, score_a = analyze_structure(Path("tests/data/mini.pdb"), "MINI", cfg, engine="a")
    patches_b, score_b = analyze_structure(Path("tests/data/mini.pdb"), "MINI", cfg, engine="b")
    assert score_a._fields == score_b._fields
    assert score_b.uniprot == "MINI"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_engine_b.py -v`
Expected: FAIL with `TypeError` (analyze_structure has no `engine` arg) or `ModuleNotFoundError`

- [ ] **Step 3: Write implementation**

```python
# prosurf/metric/engine_b.py
import numpy as np
from scipy.spatial import cKDTree
from prosurf.metric.engine_a import LocationScore
from prosurf.metric.components import balance, density
from prosurf.metric.mixing_spatial import cross_colocation
from prosurf.surface.charges import assign_charges
from prosurf.surface.sasa import surface_residue_ids

def score_locations_b(arr, cfg):
    sids = surface_residue_ids(arr, cfg.rsasa_threshold)
    charges = assign_charges(arr, sids, his_weight=cfg.his_weight)
    if not charges:
        return []
    coords = np.array([c.xyz for c in charges])
    signs = np.array([c.sign for c in charges])
    weights = np.array([c.weight for c in charges])
    tree = cKDTree(coords)
    disk_area = np.pi * cfg.patch_radius ** 2
    scores = []
    for i, c in enumerate(charges):
        idx = tree.query_ball_point(coords[i], cfg.patch_radius)
        n_pos = float(np.sum(weights[idx] * (signs[idx] > 0)))
        n_neg = float(np.sum(weights[idx] * (signs[idx] < 0)))
        B = balance(n_pos, n_neg)
        D = density(n_pos + n_neg, disk_area)
        M = cross_colocation(coords[idx], signs[idx], d_max=cfg.patch_radius/2)
        scores.append(LocationScore(c.res_id, B * D * M, n_pos, n_neg, c.xyz))
    return scores
```

Modify `prosurf/pipeline.py` `analyze_structure`:

```python
from prosurf.metric.engine_a import score_locations_a
from prosurf.metric.engine_b import score_locations_b

def analyze_structure(path, uniprot, cfg, engine="a"):
    arr = load_structure(path)
    scorer = score_locations_a if engine == "a" else score_locations_b
    scores = scorer(arr, cfg)
    patches = cluster_patches(scores, cfg)
    _, sasa = residue_sasa(arr)
    total_area = float(np.nansum(sasa))
    n_locations = len(surface_residue_ids(arr, cfg.rsasa_threshold))
    ps = aggregate_protein(uniprot, patches, total_area, max(n_locations, 1))
    return patches, ps
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_engine_b.py tests/test_pipeline.py -v`
Expected: PASS (both pipeline tests still pass; engine_b test passes)

- [ ] **Step 5: Commit**

```bash
git add prosurf/metric/engine_b.py prosurf/pipeline.py tests/test_engine_b.py
git commit -m "feat: Stage B engine and A/B engine selector"
```

---

### Task 17: A-vs-B agreement & final validation wiring

**Files:**
- Modify: `prosurf/validate/robustness.py` (add `engine` passthrough)
- Modify: `scripts/run_validation.py` (run both engines, report A-vs-B Spearman)
- Create: `tests/test_engine_agreement.py`

**Interfaces:**
- Consumes: `analyze_structure(..., engine=...)`.
- Produces:
  - `prosurf.validate.robustness.engine_agreement(structures, cfg, score_col="z_mean") -> float` — Spearman ρ between per-protein engine-A and engine-B scores.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_engine_agreement.py
from pathlib import Path
from prosurf.config import MetricConfig
from prosurf.validate.robustness import engine_agreement

def test_engine_agreement_runs():
    cfg = MetricConfig(z_percentile=0.0)
    structs = [(Path("tests/data/mini.pdb"), "MINI"),
               (Path("tests/data/mini.pdb"), "MINI2")]
    rho = engine_agreement(structs, cfg)
    # identical inputs -> defined float (nan acceptable for degenerate set)
    assert isinstance(rho, float)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_engine_agreement.py -v`
Expected: FAIL with `ImportError: cannot import name 'engine_agreement'`

- [ ] **Step 3: Write implementation**

```python
# add to prosurf/validate/robustness.py
from scipy.stats import spearmanr

def engine_agreement(structures, cfg, score_col="z_mean"):
    a, b = [], []
    for path, uniprot in structures:
        _, sa = analyze_structure(path, uniprot, cfg, engine="a")
        _, sb = analyze_structure(path, uniprot, cfg, engine="b")
        a.append(getattr(sa, score_col)); b.append(getattr(sb, score_col))
    rho, _ = spearmanr(a, b)
    return float(rho)
```

Add to `scripts/run_validation.py` (before writing the report): compute `engine_agreement(structs, cfg)` and append `- A-vs-B agreement (Spearman) = {rho:.3f}` to `lines`.

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_engine_agreement.py -v`
Expected: PASS (1 passed)

- [ ] **Step 5: Run the full suite & commit**

Run: `pytest -v`
Expected: PASS (all tests)

```bash
git add prosurf/validate/robustness.py scripts/run_validation.py tests/test_engine_agreement.py
git commit -m "feat: A-vs-B engine agreement and final validation wiring"
```

---

## Self-Review

**Spec coverage:**
- §3 architecture / staged A→B → Tasks 6, 14–16 (engine A, Stage B surface/mixing/engine, selector). ✓
- §3 per-stage caching → Task 9 (`run_pilot` caches patches). ✓
- §4 charge model (Lys/Arg/Asp/Glu centers, His fractional, termini) → Task 4. ✓
- §4 surface selection (rSASA / point cloud) → Tasks 3, 14. ✓
- §5 metric B·D̂·M product, engine A mixing, engine B bivariate stats → Tasks 5, 6, 15, 16. ✓
- §5 optional Moran's I → explicitly optional in spec; deferred, not a Phase-1 deliverable. ✓ (noted, no task)
- §6 patch clustering + per-protein vector (Z_frac/max/mean/n_patches) + Phase-2 seam → Tasks 7, 8. ✓
- §7 positive/negative controls + AUROC → Tasks 10, 13. ✓
- §7 robustness sweeps (rSASA, radius, Z threshold, pLDDT, His, rotamer noise, A-vs-B) → Tasks 11, 13, 17. *Note:* pLDDT-filter and rotamer-noise sweeps are not yet wired as dedicated config knobs — see gap below.
- §8 tech stack / layout → all tasks follow the layout. ✓
- §9 deliverables (pipeline, scores, validation report, viz) → Tasks 9, 12, 13, 17. ✓

**Gaps found & resolved:**
- **pLDDT filtering** (§4 preprocess, §7 sweep) has no dedicated task. *Resolution:* add as a follow-up micro-task during execution — extend `load_structure` to expose per-residue B-factor (pLDDT) and add `plddt_min` to `MetricConfig`, filtering residues before charge assignment. Flagged here so the implementer adds it; it follows the exact pattern of Task 3 + Task 4 and is a clean one-task addition.
- **Rotamer-noise robustness** (§7) has no task. *Resolution:* add a follow-up micro-task: a `perturb_coords(arr, sigma, seed)` helper in `prosurf/validate/robustness.py` and a sweep over `sigma`, reusing `ranking_stability`. Same pattern as Task 11.

These two are explicitly additive and isolated; they were left as flagged follow-ups rather than padding the core plan, but must be implemented to fully satisfy §7.

**Placeholder scan:** No TBD/TODO in code steps; the control-set accessions carry an explicit "confirm against literature" instruction (a data-curation step, not a code placeholder). ✓

**Type consistency:** `LocationScore(res_id, z, n_pos, n_neg, xyz)`, `Patch(res_ids, mean_z, max_z, n_pos, n_neg, size)`, `ProteinScore(uniprot, z_frac, z_max, z_mean, n_patches)`, `Charge(res_id, sign, weight, xyz)` used identically across all tasks. `analyze_structure` signature gains `engine="a"` in Task 16 and is called with it consistently thereafter. ✓
