# Cross-Species OGT Analysis — Design Spec

**Status:** Approved for planning
**Date:** 2026-06-21
**Author:** ProSurf project
**Scope:** Sub-project 1 of 2 in the cross-species expansion (see *Out of Scope* for sub-project 2).

---

## Goal

Test whether zwitterionic surface patterning (the ProSurf score **Z**) increases with
organismal **optimal growth temperature (OGT)** across the Meltome Atlas tree-of-life —
and crucially, whether it does so **over and above bulk charge composition**. The latter
is what distinguishes a genuine finding from a rediscovery of the well-known thermophile
amino-acid bias (more charged residues / ion pairs; Szilágyi & Závodszky 2000).

## Architecture

Four components, executed in order, each reusing the existing ProSurf pipeline as much as
possible:

1. **Data ingestion** — build a tidy per-protein table of `(uniprot, organism, ogt, tm)`
   from the Meltome Atlas supplementary workbooks.
2. **Scoring** — run the existing `analyze_structure` pipeline plus new composition
   features over every protein with an AlphaFold structure and ≥ 150 aa, across all species.
3. **Two-level correlation analysis** — organism-level (mean Z vs OGT) and protein-level
   (within-species Z↔Tm replication).
4. **Composition control** — partial correlation / nested regression showing Z predicts
   OGT and Tm beyond charged-residue fraction and net charge.

## Tech Stack

Python ≥ 3.11; existing deps (`biotite`, `numpy`, `scipy`, `pandas`, `requests`,
`openpyxl`, `matplotlib`). No new dependencies.

## Global Constraints

- **Reuse, don't fork** the scoring pipeline: `analyze_structure(path, uniprot, cfg, engine="a")`
  and `MetricConfig` defaults are used unchanged.
- **Consistency with the human run:** ≥ 150 aa length filter; per-organism Tm = median across
  that organism's primary (non-reproducibility) datasets; canonical UniProt ID parsing
  identical to the human run (`pid.split('_')[0].split('-')[0]`, validated against the
  canonical-accession regex).
- **Resumable:** the scoring run must checkpoint to CSV and resume after interruption,
  following the `scripts/run_all_9160.py` pattern.
- **Structures:** AlphaFold DB via the existing `fetch_af2` EBI-API path. Organisms whose
  AFDB coverage falls below 25 % of fitted proteins are reported and excluded from the
  protein-level analysis but retained in a "coverage" appendix.

---

## Background & Novelty Guard

It is already established that thermophilic proteins are charge-enriched. ProSurf's distinct
claim is about **spatial mixing of opposite charges** (the B·M part of Z = B·D̂·M), not raw
charge count. Therefore the headline statistic is **not** "Z rises with OGT" alone — it is
"**Z rises with OGT after controlling for bulk charge composition.**" Component D is the
core deliverable; components A–C set it up.

---

## Data Sources

| Source | File (cached) | Use |
|--------|---------------|-----|
| Meltome Atlas S2 (per-protein melting data, 77 sheets `ma_0001`…`ma_0077`) | `/tmp/meltome_s4.xlsx` → copy to `data/meltome_s4.xlsx` | per-protein Tm per dataset |
| Meltome Atlas S1 (dataset metadata: Dataset ID → Organism → OGT) | `/tmp/meltome_supp.xlsx` → copy to `data/meltome_s1.xlsx` | organism + OGT mapping |
| AlphaFold DB v6 | EBI prediction API | structures |

Cross-species roster (primary datasets, from the S1 *Meltome data set* sheet):

| Organism | OGT °C | Class | AFDB probe |
|----------|-------:|-------|-----------|
| Oleispira antarctica | 15 | psychrophile | ✅ 3/3 |
| Caenorhabditis elegans | 20 | mesophile | (model org) |
| Arabidopsis thaliana | 25 | mesophile | (model org) |
| Drosophila melanogaster | 28 | mesophile | (model org) |
| Danio rerio | 28 | mesophile | (model org) |
| Bacillus subtilis | 30 | mesophile | (model org) |
| Saccharomyces cerevisiae | 30 | mesophile | (model org) |
| Escherichia coli | 37 | mesophile | (model org) |
| Mus musculus | 37 | mesophile | (model org) |
| Homo sapiens | 37 | mesophile | ✅ (done, 7,505) |
| Geobacillus stearothermophilus | 55 | thermophile | ⚠️ 0/3 — likely excluded |
| Picrophilus torridus | 60 | thermophile (acidophile) | ✅ 3/3 |
| Thermus thermophilus | 70 | thermophile | ✅ 3/3 |

The human scores from `data/all_9160_scores.csv` are reused directly (re-scored only if the
composition features need backfilling — see component B).

---

## Component A — Multi-species data ingestion

**New file:** `prosurf/io/meltome.py`

**Interface (produces):**
```python
def species_datasets(s1_path: Path) -> dict[str, dict]:
    """Return {organism: {"ogt": float, "dataset_ids": [str], "class": str}}
    from the S1 'Meltome data set' sheet (primary datasets only — excludes the
    'Workflow reproducibility' and inter-lab sheets)."""

def organism_tms(s4_path: Path, dataset_ids: list[str]) -> dict[str, float]:
    """Return {canonical_uniprot: median_tm} across the given dataset sheets,
    using the same ID parsing as the human run. Proteins appearing in multiple
    of the organism's datasets take the median across them."""
```

**Script:** `scripts/build_species_table.py` →
writes `data/species_meltome_tms.csv` with columns
`uniprot, organism, ogt, organism_class, tm, n_datasets`.

## Component B — Multi-species scoring

**New file:** `prosurf/metric/composition.py`

**Interface (produces):**
```python
def composition_features(arr, cfg) -> dict:
    """Bulk charge-composition features for one structure, independent of spatial
    mixing. Returns:
      charged_frac          # (#D,E,K,R residues) / (#residues), whole protein
      surface_charged_frac  # charged surface residues / surface residues
      net_charge_per_res    # |n_pos - n_neg| / n_residues   (surface charges)
      surface_charge_density# surface charges / total surface SASA (Å^-2)
    """
```
These deliberately capture *amount* and *imbalance* of charge but NOT mixing, so they serve
as the controls that Z must beat.

**Script:** `scripts/run_multispecies.py` — generalizes `run_all_9160.py` over a species
list. For each protein: `fetch_af2` → `analyze_structure` (z_mean, z_max, z_frac, n_patches)
→ `composition_features`. Resumable checkpoint to `data/multispecies_scores.csv` with columns
`uniprot, organism, ogt, tm, z_mean, z_max, z_frac, n_patches, charged_frac,
surface_charged_frac, net_charge_per_res, surface_charge_density`. Human rows are seeded from
the existing run and backfilled with composition features in a one-time pass.

## Component C — Two-level correlation analysis

**Script:** `scripts/analyze_ogt.py`

- **Organism level:** aggregate mean and median z_mean per organism; Spearman ρ(mean Z, OGT)
  across organisms (n ≈ 12); ordinary regression with 95 % CI. Reported with the caveat that
  organism-level n is small and phylogenetically non-independent.
- **Protein level (within species):** for each organism with adequate coverage, Spearman
  ρ(z_mean, Tm) — tests whether the human ρ ≈ +0.16 generalizes. Reported as a forest-style
  table of per-organism ρ with CIs.

## Component D — Composition control (core deliverable)

Within `scripts/analyze_ogt.py`:

- **Partial correlation:** Spearman partial ρ(z_mean, OGT | charged_frac, net_charge_per_res)
  pooled across all proteins, and the analogous partial ρ(z_mean, Tm | …) within species.
- **Nested regression:** compare a composition-only model
  `Tm ~ charged_frac + net_charge_per_res + surface_charge_density`
  against `+ z_mean`; report ΔR² and the F-test / likelihood-ratio p-value for the added term,
  plus the standardized z_mean coefficient. The claim is supported iff z_mean's added
  contribution is positive and significant.

**Outputs:** `data/cross_species_report.md` + `.pdf` (reusing the WeasyPrint figure-embedding
approach already in the repo), with: OGT-trend scatter, per-organism within-species ρ forest
plot, and a bar chart of ΔR² (composition-only vs +Z).

---

## Data Flow

```
metadata workbook (S1, data/meltome_s1.xlsx) ──┐
              ├─> build_species_table.py ─> data/species_meltome_tms.csv
melting-data workbook (S2, data/meltome_s4.xlsx) ──┘                 │
                                                    ▼
                              run_multispecies.py (fetch_af2 + analyze_structure
                                                   + composition_features, resumable)
                                                    │
                                                    ▼
                                    data/multispecies_scores.csv
                                                    │
                                                    ▼
                              analyze_ogt.py ─> stats + figures
                                                 ─> data/cross_species_report.{md,pdf}
```

## Error Handling

- **Missing AFDB structure / fetch failure:** record UniProt ID in
  `data/multispecies_failed.txt`; skip; never abort the run (mirrors `run_all_9160.py`).
- **Sub-threshold organism coverage (< 25 % of fitted proteins have structures):** keep in
  organism-level OGT aggregate only if ≥ 30 proteins scored; otherwise drop and note in the
  coverage appendix. Expected casualty: *Geobacillus*.
- **Degenerate structures (no surface charges → z_mean = 0):** retained as legitimate zeros
  (as in the human run); composition features still computed.
- **Workbook parse drift (unexpected sheet/column):** `prosurf/io/meltome.py` validates the
  expected header row and raises a clear error rather than silently mis-parsing.

## Testing

- `tests/test_meltome.py` — `species_datasets` returns the 13 organisms with correct OGTs
  from a trimmed fixture workbook; `organism_tms` medians a protein appearing in two datasets.
- `tests/test_composition.py` — `composition_features` on `tests/data/mini.pdb`: charged_frac
  in [0,1]; an all-Lys synthetic gives net_charge_per_res > 0 and known charged_frac; an
  empty/charge-free structure returns zeros without error.
- `tests/test_analyze_ogt.py` — partial-correlation and nested-regression helpers on
  synthetic data where the answer is known (e.g. Z constructed to correlate with OGT only
  through a confound returns ~0 partial ρ; Z with independent signal returns > 0).

## Success Criteria

1. `data/multispecies_scores.csv` covers ≥ 10 organisms spanning OGT 15–70 °C.
2. Organism-level Spearman ρ(mean Z, OGT) reported with CI.
3. Within-species Z↔Tm ρ reported for every adequately covered organism.
4. **Core:** partial correlation and nested-regression results quantifying whether Z predicts
   OGT/Tm beyond bulk composition — with an honest verdict either way.
5. Reproducible report (`.md` + `.pdf`) committed, matching the style of the existing
   full-proteome report.

## Out of Scope (Sub-project 2)

Ortholog-pair analysis: matched thermophile↔mesophile orthologs (e.g. *Thermus* vs *E. coli*,
*Picrophilus* vs yeast), ΔZ vs ΔTm within pairs. Requires new ortholog-mapping infrastructure
(OrthoDB / eggNOG / reciprocal-best-hit) and gets its own spec after sub-project 1 ships.
