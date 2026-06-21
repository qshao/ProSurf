# Cross-Species OGT Analysis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Test whether the ProSurf zwitterionic surface score (Z) rises with organismal optimal growth temperature (OGT) across the Meltome Atlas tree-of-life, controlling for bulk charge composition.

**Architecture:** Three new library modules (Meltome workbook ingestion, bulk-composition features, cross-species statistics) feed three orchestration scripts (build species Tm table → resumable multi-species scoring → analysis + report). The existing `analyze_structure` pipeline and `MetricConfig` are reused unchanged.

**Tech Stack:** Python ≥ 3.11; existing deps only — `biotite`, `numpy`, `scipy`, `pandas`, `requests`, `openpyxl`, `matplotlib`, `weasyprint`. No new dependencies.

## Global Constraints

- Reuse `analyze_structure(path, uniprot, cfg, engine="a")` and `MetricConfig()` defaults unchanged.
- Length filter: skip proteins < 150 aa.
- Per-organism Tm = median across that organism's primary (non-reproducibility) datasets.
- Canonical UniProt ID parsing: `pid.split('_')[0].split('-')[0]`, validated against regex `^([A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]|[OPQ][0-9][A-Z0-9]{3}[0-9])$`.
- Scoring run must checkpoint to CSV and resume after interruption.
- Structures via `fetch_af2` (EBI API). Fetch/parse failures are logged and skipped, never abort the run.
- Headline statistic is Z's contribution **beyond** bulk composition (partial correlation + nested regression), not the raw OGT trend alone.

---

### Task 1: Meltome workbook ingestion module

**Files:**
- Create: `prosurf/io/meltome.py`
- Test: `tests/test_meltome.py`
- Test fixtures: created in the test from openpyxl (no binary fixtures committed)

**Interfaces:**
- Consumes: nothing (entry layer).
- Produces:
  - `canonical_uniprot(protein_id: str) -> str | None`
  - `species_datasets(s1_path: Path) -> dict[str, dict]` — `{organism: {"ogt": float, "dataset_ids": list[str], "class": str}}`
  - `organism_tms(s4_path: Path, dataset_ids: list[str]) -> dict[str, float]` — `{canonical_uniprot: median_tm}`

- [ ] **Step 1: Write the failing test**

Create `tests/test_meltome.py`:

```python
from pathlib import Path
import openpyxl
from prosurf.io.meltome import canonical_uniprot, species_datasets, organism_tms


def test_canonical_uniprot_parses_and_filters():
    assert canonical_uniprot("P61626_LYZ") == "P61626"
    assert canonical_uniprot("Q9UNM6-2_PSMD13") == "Q9UNM6"
    assert canonical_uniprot("O50146_lysY") == "O50146"
    assert canonical_uniprot("notanid_xyz") is None


def _make_s1(path):
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Meltome data set"
    ws.append(["Dataset ID", "Dataset source", "Organism", "Strain/tissue", "OGT [°C]"])
    ws.append(["ma:0008", "TUM", "Thermus thermophilus", "HB27", 70])
    ws.append(["ma:0009", "TUM", "Thermus thermophilus", "HB27", 70])
    ws.append(["ma:0003", "TUM", "Escherichia coli", "K12", 37])
    ws.append(["ma:0001", "TUM", "Oleispira antarctica", "RB-8", 15])
    wb.save(path)


def _make_s4(path):
    wb = openpyxl.Workbook()
    first = True
    for sheet in ["ma_0008", "ma_0009"]:
        ws = wb.create_sheet(sheet) if not first else wb.active
        if first:
            ws.title = sheet
            first = False
        ws.append(["Protein ID", "Melting point [°C]", "class"])
    # ma_0008: O50147 melts at 80, P0 non-melter (None)
    wb["ma_0008"].append(["O50147_lysZ", 80.0, "medium"])
    wb["ma_0008"].append(["O66271_fumC", None, "non-melter"])
    # ma_0009: O50147 again at 82 -> median across datasets = 81
    wb["ma_0009"].append(["O50147_lysZ", 82.0, "medium"])
    wb.save(path)


def test_species_datasets(tmp_path):
    p = tmp_path / "s1.xlsx"
    _make_s1(p)
    spec = species_datasets(p)
    assert spec["Thermus thermophilus"]["ogt"] == 70.0
    assert sorted(spec["Thermus thermophilus"]["dataset_ids"]) == ["ma_0008", "ma_0009"]
    assert spec["Thermus thermophilus"]["class"] == "thermophile"
    assert spec["Escherichia coli"]["class"] == "mesophile"
    assert spec["Oleispira antarctica"]["class"] == "psychrophile"


def test_organism_tms_medians_across_datasets(tmp_path):
    p = tmp_path / "s4.xlsx"
    _make_s4(p)
    tms = organism_tms(p, ["ma_0008", "ma_0009"])
    assert tms["O50147"] == 81.0      # median of 80 and 82
    assert "O66271" not in tms        # non-melter (None Tm) dropped
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_meltome.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'prosurf.io.meltome'`

- [ ] **Step 3: Write minimal implementation**

Create `prosurf/io/meltome.py`:

```python
import re
import statistics
from pathlib import Path

import openpyxl

_CANONICAL_RE = re.compile(
    r"^([A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]|[OPQ][0-9][A-Z0-9]{3}[0-9])$"
)


def _organism_class(ogt: float) -> str:
    if ogt < 20:
        return "psychrophile"
    if ogt <= 45:
        return "mesophile"
    return "thermophile"


def canonical_uniprot(protein_id: str) -> str | None:
    """Extract canonical UniProt accession from a Meltome protein ID such as
    'P61626_LYZ' or 'Q9UNM6-2_PSMD13'. Returns None if not canonical."""
    uid = protein_id.split("_")[0].split("-")[0]
    return uid if _CANONICAL_RE.match(uid) else None


def species_datasets(s1_path: Path) -> dict:
    """Parse the S1 'Meltome data set' sheet into
    {organism: {'ogt': float, 'dataset_ids': [str], 'class': str}}.
    Only the primary 'Meltome data set' sheet is read; reproducibility and
    inter-lab sheets are ignored."""
    wb = openpyxl.load_workbook(s1_path, read_only=True)
    ws = wb["Meltome data set"]
    rows = list(ws.iter_rows(values_only=True))
    header = rows[0]
    if header[0] != "Dataset ID" or "OGT [°C]" not in header:
        raise ValueError(f"Unexpected S1 header: {header[:5]}")
    org_idx = header.index("Organism")
    ogt_idx = header.index("OGT [°C]")
    out: dict = {}
    for r in rows[1:]:
        if not r or not r[0]:
            continue
        did = str(r[0]).replace("ma:", "ma_")
        org = r[org_idx]
        ogt = r[ogt_idx]
        if org is None or not isinstance(ogt, (int, float)):
            continue
        d = out.setdefault(
            org,
            {"ogt": float(ogt), "dataset_ids": [], "class": _organism_class(float(ogt))},
        )
        d["dataset_ids"].append(did)
    return out


def organism_tms(s4_path: Path, dataset_ids: list) -> dict:
    """Median Tm per canonical UniProt across the given S4 dataset sheets.
    Sheet layout: column 0 = Protein ID, column 1 = Melting point [°C]
    (None for non-melters, which are dropped)."""
    wb = openpyxl.load_workbook(s4_path, read_only=True)
    per_protein: dict = {}
    for did in dataset_ids:
        if did not in wb.sheetnames:
            continue
        ws = wb[did]
        for i, row in enumerate(ws.iter_rows(values_only=True)):
            if i == 0 or not row or not row[0]:
                continue
            uid = canonical_uniprot(str(row[0]))
            tm = row[1] if len(row) > 1 else None
            if uid is None or not isinstance(tm, (int, float)):
                continue
            per_protein.setdefault(uid, []).append(float(tm))
    return {uid: statistics.median(v) for uid, v in per_protein.items()}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_meltome.py -v`
Expected: PASS (4 tests)

- [ ] **Step 5: Commit**

```bash
git add prosurf/io/meltome.py tests/test_meltome.py
git commit -m "feat: Meltome workbook ingestion (species, OGT, per-organism Tm)"
```

---

### Task 2: Species Tm table build script

**Files:**
- Create: `scripts/build_species_table.py`
- Copy workbooks into repo: `data/meltome_s1.xlsx`, `data/meltome_s4.xlsx`
- Output: `data/species_meltome_tms.csv`

**Interfaces:**
- Consumes: `prosurf.io.meltome.species_datasets`, `prosurf.io.meltome.organism_tms`.
- Produces: `data/species_meltome_tms.csv` with header `uniprot,organism,ogt,organism_class,tm,n_datasets`.

- [ ] **Step 1: Copy the cached workbooks into the repo**

```bash
cp /tmp/meltome_supp.xlsx data/meltome_s1.xlsx
cp /tmp/meltome_s4.xlsx data/meltome_s4.xlsx
ls -l data/meltome_s1.xlsx data/meltome_s4.xlsx
```
Expected: both files listed (≈ 32 KB and ≈ 14 MB).

- [ ] **Step 2: Write the build script**

Create `scripts/build_species_table.py`:

```python
"""Build data/species_meltome_tms.csv from the Meltome Atlas workbooks.

One row per (protein, organism): uniprot, organism, ogt, organism_class,
tm (median across that organism's primary datasets), n_datasets.
Human is included here too; the scoring step reuses existing human scores.
"""
import csv
from pathlib import Path

from prosurf.io.meltome import species_datasets, organism_tms

ROOT = Path(__file__).resolve().parent.parent
S1 = ROOT / "data" / "meltome_s1.xlsx"
S4 = ROOT / "data" / "meltome_s4.xlsx"
OUT = ROOT / "data" / "species_meltome_tms.csv"


def main():
    spec = species_datasets(S1)
    print(f"Found {len(spec)} organisms in S1")
    rows = []
    for organism, info in sorted(spec.items(), key=lambda kv: kv[1]["ogt"]):
        tms = organism_tms(S4, info["dataset_ids"])
        print(f"  {organism:38s} OGT={info['ogt']:>4}  "
              f"datasets={len(info['dataset_ids'])}  proteins={len(tms)}")
        for uid, tm in tms.items():
            rows.append({
                "uniprot": uid,
                "organism": organism,
                "ogt": info["ogt"],
                "organism_class": info["class"],
                "tm": round(tm, 3),
                "n_datasets": len(info["dataset_ids"]),
            })
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=["uniprot", "organism", "ogt",
                           "organism_class", "tm", "n_datasets"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {len(rows)} rows to {OUT}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Run the script**

Run: `python3 scripts/build_species_table.py`
Expected: prints ~13 organisms ordered by OGT (Oleispira 15 … Thermus 70), then `Wrote N rows to data/species_meltome_tms.csv` with N in the tens of thousands.

- [ ] **Step 4: Verify the output shape**

Run:
```bash
python3 -c "
import pandas as pd
df = pd.read_csv('data/species_meltome_tms.csv')
print('rows', len(df))
print('organisms', df.organism.nunique())
print('ogt range', df.ogt.min(), df.ogt.max())
assert {'uniprot','organism','ogt','organism_class','tm','n_datasets'} <= set(df.columns)
assert df.ogt.min() <= 15 and df.ogt.max() >= 70
print('OK')
"
```
Expected: ends with `OK`.

- [ ] **Step 5: Commit**

```bash
git add scripts/build_species_table.py data/meltome_s1.xlsx data/meltome_s4.xlsx data/species_meltome_tms.csv
git commit -m "feat: build per-organism Tm table from Meltome workbooks"
```

---

### Task 3: Bulk-composition features

**Files:**
- Create: `prosurf/metric/composition.py`
- Test: `tests/test_composition.py`

**Interfaces:**
- Consumes: `prosurf.surface.sasa.residue_sasa`, `prosurf.surface.sasa.surface_residue_ids`, `prosurf.surface.charges.assign_charges`, `MetricConfig`.
- Produces: `composition_features(arr, cfg) -> dict` with keys `charged_frac`, `surface_charged_frac`, `net_charge_per_res`, `surface_charge_density` (all float).

- [ ] **Step 1: Write the failing test**

Create `tests/test_composition.py`:

```python
from pathlib import Path
from prosurf.config import MetricConfig
from prosurf.io.parse import load_structure
from prosurf.metric.composition import composition_features


def test_composition_on_mini():
    # mini.pdb = LYS, ASP, ALA -> 2 of 3 residues charged
    arr = load_structure(Path("tests/data/mini.pdb"))
    cfg = MetricConfig()
    feats = composition_features(arr, cfg)
    assert abs(feats["charged_frac"] - 2 / 3) < 1e-6
    for k in ("surface_charged_frac", "net_charge_per_res", "surface_charge_density"):
        assert feats[k] >= 0.0
    assert 0.0 <= feats["surface_charged_frac"] <= 1.0


def test_composition_keys_are_floats():
    arr = load_structure(Path("tests/data/mini.pdb"))
    feats = composition_features(arr, MetricConfig())
    assert set(feats) == {"charged_frac", "surface_charged_frac",
                          "net_charge_per_res", "surface_charge_density"}
    assert all(isinstance(v, float) for v in feats.values())
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_composition.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'prosurf.metric.composition'`

- [ ] **Step 3: Write minimal implementation**

Create `prosurf/metric/composition.py`:

```python
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_composition.py -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Commit**

```bash
git add prosurf/metric/composition.py tests/test_composition.py
git commit -m "feat: bulk charge-composition features (control covariates for Z)"
```

---

### Task 4: Resumable multi-species scoring script

**Files:**
- Create: `scripts/run_multispecies.py`
- Output: `data/multispecies_scores.csv`, `data/multispecies_failed.txt`

**Interfaces:**
- Consumes: `data/species_meltome_tms.csv` (Task 2); `prosurf.io.fetch.fetch_af2`; `prosurf.pipeline.analyze_structure`; `prosurf.metric.composition.composition_features`; `prosurf.config.MetricConfig`, `PathsConfig`; `prosurf.io.parse.load_structure`.
- Produces: `data/multispecies_scores.csv` with header `uniprot,organism,ogt,tm,z_mean,z_max,z_frac,n_patches,charged_frac,surface_charged_frac,net_charge_per_res,surface_charge_density`.

- [ ] **Step 1: Write the script**

Create `scripts/run_multispecies.py`:

```python
"""Score every protein in data/species_meltome_tms.csv across all organisms.

Resumable: appends to data/multispecies_scores.csv and skips proteins already
present there or in data/multispecies_failed.txt. Length filter < 150 aa.

Flags:
  --organisms "A,B"   restrict to these organisms (substring match)
  --limit N           cap proteins scored per organism (smoke testing)
  --min-len N         length floor (default 150)
"""
import argparse
import csv
import time
from pathlib import Path

import biotite.structure as struc

from prosurf.config import MetricConfig, PathsConfig
from prosurf.io.fetch import fetch_af2
from prosurf.io.parse import load_structure
from prosurf.metric.composition import composition_features
from prosurf.pipeline import analyze_structure

ROOT = Path(__file__).resolve().parent.parent
IN_CSV = ROOT / "data" / "species_meltome_tms.csv"
OUT_CSV = ROOT / "data" / "multispecies_scores.csv"
FAILED = ROOT / "data" / "multispecies_failed.txt"

FIELDS = ["uniprot", "organism", "ogt", "tm", "z_mean", "z_max", "z_frac",
          "n_patches", "charged_frac", "surface_charged_frac",
          "net_charge_per_res", "surface_charge_density"]


def load_done():
    done = set()
    if OUT_CSV.exists():
        with open(OUT_CSV) as f:
            for r in csv.DictReader(f):
                done.add((r["uniprot"], r["organism"]))
    return done


def load_failed():
    if FAILED.exists():
        return set(FAILED.read_text().splitlines())
    return set()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organisms", default="")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--min-len", type=int, default=150)
    args = ap.parse_args()

    cfg = MetricConfig()
    paths = PathsConfig()

    with open(IN_CSV) as f:
        targets = list(csv.DictReader(f))
    if args.organisms:
        wanted = [s.strip().lower() for s in args.organisms.split(",")]
        targets = [t for t in targets
                   if any(w in t["organism"].lower() for w in wanted)]

    done = load_done()
    failed = load_failed()

    # Seed existing human scores from the full-proteome run (reuse, don't refetch).
    human_seed = ROOT / "data" / "all_9160_scores.csv"
    seeded = {}
    if human_seed.exists():
        with open(human_seed) as f:
            for r in csv.DictReader(f):
                seeded[r["uniprot"]] = r

    if not OUT_CSV.exists():
        with open(OUT_CSV, "w", newline="") as f:
            csv.DictWriter(f, fieldnames=FIELDS).writeheader()

    # Per-organism caps for --limit
    per_org_count: dict = {}
    t0 = time.time()
    n_new = 0

    for t in targets:
        uid, organism = t["uniprot"], t["organism"]
        if (uid, organism) in done or uid in failed:
            continue
        if args.limit:
            c = per_org_count.get(organism, 0)
            if c >= args.limit:
                continue

        try:
            pdb = fetch_af2(uid, paths.data_dir)
            arr = load_structure(pdb)
            n_res = len(struc.get_residues(arr)[0])
            if n_res < args.min_len:
                with open(FAILED, "a") as ff:
                    ff.write(uid + "\n")
                continue
            _, ps = analyze_structure(pdb, uid, cfg)
            comp = composition_features(arr, cfg)
            row = {
                "uniprot": uid, "organism": organism, "ogt": t["ogt"], "tm": t["tm"],
                "z_mean": ps.z_mean, "z_max": ps.z_max, "z_frac": ps.z_frac,
                "n_patches": ps.n_patches, **comp,
            }
            with open(OUT_CSV, "a", newline="") as f:
                csv.DictWriter(f, fieldnames=FIELDS).writerow(row)
            n_new += 1
            per_org_count[organism] = per_org_count.get(organism, 0) + 1
            if n_new % 25 == 0:
                rate = n_new / (time.time() - t0)
                print(f"  {n_new} scored  ({organism} {uid} "
                      f"z_mean={ps.z_mean:.3f})  {rate:.1f}/s", flush=True)
        except Exception as e:
            with open(FAILED, "a") as ff:
                ff.write(uid + "\n")
            print(f"  FAIL {organism} {uid}: {e}", flush=True)

    print(f"Done. {n_new} new proteins scored -> {OUT_CSV}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Smoke-test on one small organism**

Run: `python3 scripts/run_multispecies.py --organisms "Thermus" --limit 10`
Expected: scores up to 10 *Thermus thermophilus* proteins, prints `Done. N new proteins scored` with N ≤ 10, no traceback.

- [ ] **Step 3: Verify output columns and resume behaviour**

Run:
```bash
python3 -c "
import pandas as pd
df = pd.read_csv('data/multispecies_scores.csv')
expected = ['uniprot','organism','ogt','tm','z_mean','z_max','z_frac','n_patches','charged_frac','surface_charged_frac','net_charge_per_res','surface_charge_density']
assert list(df.columns) == expected, df.columns.tolist()
assert (df.organism == 'Thermus thermophilus').all()
print('rows', len(df), 'OK')
"
python3 scripts/run_multispecies.py --organisms "Thermus" --limit 10
```
Expected: first command prints `rows N OK`; the re-run prints `Done. 0 new proteins scored` (all already done — resume works).

- [ ] **Step 4: Commit**

```bash
git add scripts/run_multispecies.py data/multispecies_scores.csv data/multispecies_failed.txt
git commit -m "feat: resumable multi-species scoring with composition features"
```

- [ ] **Step 5: Launch the full run (long-running, background)**

```bash
tmux new-session -d -s prosurf_species \
  'cd '"$(pwd)"' && python3 scripts/run_multispecies.py 2>&1 | tee logs/run_multispecies.log'
```
Expected: returns immediately; `tail -f logs/run_multispecies.log` shows progress. Let it finish (all organisms) before Task 6. Confirm completion with `tail -1 logs/run_multispecies.log` showing `Done. N new proteins scored`.

---

### Task 5: Cross-species statistics helpers

**Files:**
- Create: `prosurf/validate/crossspecies.py`
- Test: `tests/test_crossspecies.py`

**Interfaces:**
- Consumes: `numpy`, `scipy.stats`, `pandas`.
- Produces:
  - `partial_spearman(x, y, covariates) -> tuple[float, float]` — (partial r, p-value); `covariates` is a list of 1-D array-likes.
  - `nested_regression_delta(df, target, base_features, added_feature) -> dict` with keys `r2_base`, `r2_full`, `delta_r2`, `f_stat`, `p_value`, `added_coef`.

- [ ] **Step 1: Write the failing test**

Create `tests/test_crossspecies.py`:

```python
import numpy as np
import pandas as pd
from prosurf.validate.crossspecies import partial_spearman, nested_regression_delta


def test_partial_spearman_removes_confound():
    rng = np.random.default_rng(0)
    z = rng.normal(size=400)              # confound
    x = z + 0.01 * rng.normal(size=400)   # x driven only by z
    y = z + 0.01 * rng.normal(size=400)   # y driven only by z
    # raw correlation high, partial (controlling z) near zero
    r_partial, p = partial_spearman(x, y, [z])
    assert abs(r_partial) < 0.2


def test_partial_spearman_keeps_real_signal():
    rng = np.random.default_rng(1)
    z = rng.normal(size=400)
    sig = rng.normal(size=400)
    x = z + sig
    y = z + sig                            # shared signal beyond z
    r_partial, p = partial_spearman(x, y, [z])
    assert r_partial > 0.4
    assert p < 1e-3


def test_nested_regression_delta_detects_added_signal():
    rng = np.random.default_rng(2)
    base = rng.normal(size=500)
    extra = rng.normal(size=500)
    y = 1.0 * base + 1.0 * extra + 0.1 * rng.normal(size=500)
    df = pd.DataFrame({"y": y, "base": base, "extra": extra})
    res = nested_regression_delta(df, "y", ["base"], "extra")
    assert res["delta_r2"] > 0.2
    assert res["p_value"] < 1e-6
    assert res["added_coef"] > 0


def test_nested_regression_delta_null_added_signal():
    rng = np.random.default_rng(3)
    base = rng.normal(size=500)
    noise = rng.normal(size=500)           # unrelated to y
    y = base + 0.1 * rng.normal(size=500)
    df = pd.DataFrame({"y": y, "base": base, "noise": noise})
    res = nested_regression_delta(df, "y", ["base"], "noise")
    assert res["delta_r2"] < 0.02
    assert res["p_value"] > 0.05
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pytest tests/test_crossspecies.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'prosurf.validate.crossspecies'`

- [ ] **Step 3: Write minimal implementation**

Create `prosurf/validate/crossspecies.py`:

```python
import numpy as np
from scipy import stats


def partial_spearman(x, y, covariates):
    """Partial Spearman correlation of x and y controlling for covariates.

    Rank-transforms all variables, regresses ranked x and ranked y on the
    ranked covariates (with intercept), and correlates the residuals.
    `covariates` is a list of 1-D array-likes. Returns (partial_r, p_value).
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    C = np.column_stack([np.asarray(c, float) for c in covariates])

    xr = stats.rankdata(x)
    yr = stats.rankdata(y)
    Cr = np.column_stack([stats.rankdata(C[:, j]) for j in range(C.shape[1])])
    A = np.column_stack([np.ones(len(xr)), Cr])

    bx, *_ = np.linalg.lstsq(A, xr, rcond=None)
    by, *_ = np.linalg.lstsq(A, yr, rcond=None)
    rx = xr - A @ bx
    ry = yr - A @ by
    r, p = stats.pearsonr(rx, ry)
    return float(r), float(p)


def nested_regression_delta(df, target, base_features, added_feature):
    """Compare OLS models target ~ base vs target ~ base + added_feature.

    Predictors are z-standardized. Returns dict: r2_base, r2_full, delta_r2,
    f_stat, p_value (F-test for the single added term), added_coef (standardized).
    """
    y = df[target].to_numpy(float)
    n = len(y)

    def standardize(M):
        M = np.asarray(M, float)
        return (M - M.mean(0)) / M.std(0)

    Xb = standardize(df[base_features].to_numpy(float))
    Xf = standardize(df[base_features + [added_feature]].to_numpy(float))
    Ab = np.column_stack([np.ones(n), Xb])
    Af = np.column_stack([np.ones(n), Xf])

    def fit(A):
        beta, *_ = np.linalg.lstsq(A, y, rcond=None)
        resid = y - A @ beta
        rss = float(resid @ resid)
        tss = float(((y - y.mean()) ** 2).sum())
        return beta, rss, 1.0 - rss / tss

    bb, rss_b, r2_b = fit(Ab)
    bf, rss_f, r2_f = fit(Af)

    df_num = Af.shape[1] - Ab.shape[1]
    df_den = n - Af.shape[1]
    f_stat = ((rss_b - rss_f) / df_num) / (rss_f / df_den)
    p_value = float(stats.f.sf(f_stat, df_num, df_den))

    return {
        "r2_base": r2_b, "r2_full": r2_f, "delta_r2": r2_f - r2_b,
        "f_stat": float(f_stat), "p_value": p_value, "added_coef": float(bf[-1]),
    }
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pytest tests/test_crossspecies.py -v`
Expected: PASS (4 tests)

- [ ] **Step 5: Commit**

```bash
git add prosurf/validate/crossspecies.py tests/test_crossspecies.py
git commit -m "feat: partial correlation + nested regression for cross-species control"
```

---

### Task 6: OGT analysis and report

**Files:**
- Create: `scripts/analyze_ogt.py`
- Output: `data/cross_species_report.md`, `data/cross_species_report.pdf`

**Interfaces:**
- Consumes: `data/multispecies_scores.csv` (Task 4); `prosurf.validate.crossspecies.partial_spearman`, `nested_regression_delta`; `scipy.stats.spearmanr`; `matplotlib`; `weasyprint`; `markdown`.
- Produces: the report files. No new importable API.

**Precondition:** Task 4's full run has completed (`data/multispecies_scores.csv` covers ≥ 10 organisms).

- [ ] **Step 1: Write the analysis script**

Create `scripts/analyze_ogt.py`:

```python
"""Cross-species OGT analysis: organism-level trend, within-species replication,
and the composition control. Writes data/cross_species_report.{md,pdf}."""
import base64
import io
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from prosurf.validate.crossspecies import partial_spearman, nested_regression_delta

ROOT = Path(__file__).resolve().parent.parent
SCORES = ROOT / "data" / "multispecies_scores.csv"
MD = ROOT / "data" / "cross_species_report.md"
PDF = ROOT / "data" / "cross_species_report.pdf"

COVARS = ["charged_frac", "net_charge_per_res", "surface_charge_density"]
MIN_PER_ORG = 30


def fig_b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def main():
    df = pd.read_csv(SCORES)
    df = df.dropna(subset=["z_mean", "tm", "ogt"] + COVARS)

    # Organism-level trend (drop organisms with < MIN_PER_ORG scored proteins)
    org = (df.groupby("organism")
             .agg(ogt=("ogt", "first"), mean_z=("z_mean", "mean"),
                  n=("z_mean", "size"))
             .reset_index())
    org = org[org["n"] >= MIN_PER_ORG].sort_values("ogt")
    rho_org, p_org = spearmanr(org["mean_z"], org["ogt"])

    # Within-species replication
    within = []
    for organism, g in df.groupby("organism"):
        if len(g) >= MIN_PER_ORG:
            r, p = spearmanr(g["z_mean"], g["tm"])
            within.append((organism, g["ogt"].iloc[0], len(g), r, p))
    within.sort(key=lambda t: t[1])

    # Composition control (pooled): Z vs OGT controlling for bulk composition
    pr_ogt, pp_ogt = partial_spearman(
        df["z_mean"], df["ogt"], [df[c] for c in COVARS])
    raw_ogt, _ = spearmanr(df["z_mean"], df["ogt"])

    # Nested regression: does Z add to a composition-only model of Tm?
    nested = nested_regression_delta(df, "tm", COVARS, "z_mean")

    # ---- Figures ----
    f1, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(org["ogt"], org["mean_z"], s=60, color="#2c3e50", zorder=3)
    for _, r in org.iterrows():
        ax.annotate(r["organism"].split()[0], (r["ogt"], r["mean_z"]),
                    fontsize=7, xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel("Optimal growth temperature (°C)")
    ax.set_ylabel("Mean z_mean")
    ax.set_title(f"Organism-level: mean Z vs OGT\nSpearman ρ={rho_org:+.3f}, p={p_org:.3f}")
    ax.grid(True, lw=0.4, color="#ccc")
    img_org = fig_b64(f1)

    f2, ax = plt.subplots(figsize=(7, max(3, 0.4 * len(within) + 1)))
    labels = [w[0] for w in within]
    rs = [w[3] for w in within]
    ax.barh(range(len(within)), rs,
            color=["#27ae60" if r > 0 else "#e74c3c" for r in rs])
    ax.set_yticks(range(len(within)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(0, color="#333", lw=0.8)
    ax.set_xlabel("within-species Spearman ρ(z_mean, Tm)")
    ax.set_title("Within-species replication")
    img_within = fig_b64(f2)

    f3, ax = plt.subplots(figsize=(5, 4.5))
    ax.bar(["composition\nonly", "+ z_mean"],
           [nested["r2_base"], nested["r2_full"]],
           color=["#95a5a6", "#2c3e50"])
    ax.set_ylabel("R² (model of Tm)")
    ax.set_title(f"Nested model: ΔR²={nested['delta_r2']:+.4f}\n"
                 f"added-term p={nested['p_value']:.2e}")
    img_delta = fig_b64(f3)

    # ---- Markdown ----
    within_rows = "\n".join(
        f"| {o} | {ogt:.0f} | {n} | {r:+.3f} | {p:.2e} |"
        for o, ogt, n, r, p in within)
    md = f"""# Cross-Species OGT Analysis

**N proteins:** {len(df):,} across {df.organism.nunique()} organisms (OGT {df.ogt.min():.0f}–{df.ogt.max():.0f} °C)

## 1. Organism-level trend
Spearman ρ(mean Z, OGT) = **{rho_org:+.3f}** (p = {p_org:.3f}, n = {len(org)} organisms).

![org](data:image/png;base64,{img_org})

## 2. Within-species replication

| Organism | OGT | N | ρ(z_mean, Tm) | p |
|----------|----:|--:|--------------:|---|
{within_rows}

![within](data:image/png;base64,{img_within})

## 3. Composition control (headline)
- Raw Spearman ρ(z_mean, OGT) = {raw_ogt:+.3f}
- **Partial** ρ(z_mean, OGT | {', '.join(COVARS)}) = **{pr_ogt:+.3f}** (p = {pp_ogt:.2e})
- Nested model of Tm: R²(composition) = {nested['r2_base']:.4f} → R²(+z_mean) = {nested['r2_full']:.4f}; **ΔR² = {nested['delta_r2']:+.4f}**, added-term p = {nested['p_value']:.2e}, standardized z_mean coef = {nested['added_coef']:+.3f}.

![delta](data:image/png;base64,{img_delta})

**Verdict:** Z {'retains' if pr_ogt > 0 and pp_ogt < 0.05 else 'does NOT retain'} a significant positive association with thermal stability after controlling for bulk charge composition.
"""
    MD.write_text(md)
    print(f"Wrote {MD}")
    print(f"  organism-level rho={rho_org:+.3f} p={p_org:.3f}")
    print(f"  partial rho(Z,OGT|comp)={pr_ogt:+.3f} p={pp_ogt:.2e}")
    print(f"  nested deltaR2={nested['delta_r2']:+.4f} p={nested['p_value']:.2e}")

    # ---- PDF ----
    try:
        import markdown
        from weasyprint import HTML
        html = markdown.markdown(md, extensions=["tables"])
        HTML(string=html).write_pdf(PDF)
        print(f"Wrote {PDF}")
    except Exception as e:
        print(f"PDF generation skipped: {e}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the analysis**

Run: `python3 scripts/analyze_ogt.py`
Expected: prints `Wrote data/cross_species_report.md`, the three summary lines (organism-level rho, partial rho, nested deltaR2), and `Wrote data/cross_species_report.pdf`.

- [ ] **Step 3: Sanity-check the report**

Run:
```bash
python3 -c "
import pathlib
t = pathlib.Path('data/cross_species_report.md').read_text()
assert 'Composition control' in t
assert 'Within-species replication' in t
assert 'Verdict' in t
print('report OK,', len(t), 'chars')
"
ls -l data/cross_species_report.pdf
```
Expected: `report OK, N chars` and the PDF file listed with non-zero size.

- [ ] **Step 4: Commit**

```bash
git add scripts/analyze_ogt.py data/cross_species_report.md data/cross_species_report.pdf
git commit -m "feat: cross-species OGT analysis report (trend, replication, composition control)"
```

---

## Final verification

- [ ] Run the full test suite: `pytest tests/ -q` — expected: all pass.
- [ ] Confirm the four deliverable data files exist: `data/species_meltome_tms.csv`, `data/multispecies_scores.csv`, `data/cross_species_report.md`, `data/cross_species_report.pdf`.
- [ ] Confirm `data/multispecies_scores.csv` covers ≥ 10 organisms spanning OGT 15–70 °C.
