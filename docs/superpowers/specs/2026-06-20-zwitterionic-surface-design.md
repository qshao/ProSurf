# Zwitterionic Surface Pattern Quantification — Design Spec

**Date:** 2026-06-20
**Project:** ProSurf
**Status:** Phase 1 design (algorithm development)

## 1. Motivation & Hypothesis

**Hypothesis:** Protein surfaces that present a *zwitterionic pattern* tend to
have higher stability.

A **zwitterionic surface region** is defined here as a surface area that is
**net-charge neutral but charge-rich** — i.e., it achieves local neutrality
through an intermixing of *both* positive (Lys/Arg) and negative (Asp/Glu)
residues in close proximity, possibly together with other hydrophilic or
hydrophobic residues. This is explicitly different from a region that is
neutral simply because it lacks charges, and different from a region where +
and − charges are segregated into separate sub-domains.

This spec covers **Phase 1**: develop and validate a robust algorithm to
**quantify the zwitterionic pattern** on AlphaFold2-predicted human protein
structures. Correlation with a stability signal is **Phase 2** and is
deliberately deferred; the design exposes a clean seam where the stability
signal attaches with no upstream rework.

## 2. Scope

- **Structures:** AlphaFold2-predicted human protein structures.
- **Phase 1 scale:** a curated pilot subset of ~100–1000 proteins (diverse,
  including the validation control set). Proteome-wide scaling is a later phase.
- **Primary output:** *both* (a) a per-patch map locating zwitterionic regions
  on each surface, and (b) a per-protein score vector aggregated from the map.
  The patch map is foundational; the protein scores derive from it.
- **Out of scope (Phase 1):** stability data, stability correlation, full
  proteome run.

## 3. Architecture & Data Flow

A staged, modular pipeline. Each stage reads/writes well-defined on-disk
artifacts so stages are independently testable, cacheable, and re-runnable.

```
[1] Acquire        AF2 human structures (PDB/mmCIF) + metadata
       │
[2] Preprocess     parse, select chains, pLDDT filter, compute SASA
       │
[3] Surface model  surface residues (A) / SES mesh (B); assign charges
       │
[4] Zwitterion     per-location Z score
    engine           A: composite B·D̂·M
                     B: bivariate spatial statistics
       │
[5] Patches        cluster high-Z locations → patch map per protein
       │
[6] Aggregate      per-protein score vector  (← Phase-2 stability hook)
       │
[7] Validate/viz   controls, robustness sweeps, render maps
```

**Two interchangeable backends** sit behind one interface for stages [3]+[4]:
`engine="A"` (residue / Euclidean sphere) and `engine="B"` (mesh / geodesic).
Everything downstream consumes the same **`Z-per-location`** artifact, so the
A→B transition is a backend swap, not a rewrite.

**Key design choices:**
- **Per-stage caching** keyed by input + parameters, so robustness sweeps
  re-run only the affected stage.
- **One protein = one independent unit of work**, trivially parallelizable for
  later scaling.

### Staged development plan (A → B)
- **Stage A (prototype, ~1 week target):** residue-graph + spherical patches +
  composite index. Nail mechanics, visualization, and the control set fast.
- **Stage B (production metric):** molecular-surface mesh + geodesic patches +
  bivariate spatial statistics. Becomes the primary metric.
- C (Poisson–Boltzmann electrostatic field) is **not** in scope but may be
  added later as a continuous physics-based cross-check.

## 4. Surface & Charge Model

Consistent inputs shared by both engines.

**Charge assignment (physiological pH ≈ 7.4):**
- **Positive:** Lys (point charge at NZ), Arg (at guanidinium center, CZ).
- **Negative:** Asp (carboxylate center, midpoint of OD1/OD2), Glu (midpoint of
  OE1/OE2).
- **His:** **excluded by default** (pKa ≈ 6 → ~5–10% protonated at pH 7.4),
  exposed as a tunable **fractional weight** (0 → 0.1) for robustness checks.
- **Termini:** N-terminus (+) and C-terminus (−) included as point charges at
  the terminal atoms.
- Each charge is a **point at its side-chain charge center** (not Cα), since the
  charged group projects several Å from the backbone.

**Surface selection:**
- *Engine A:* residue is "surface" if relative SASA > threshold (default 20%,
  swept in robustness).
- *Engine B:* the solvent-excluded surface (SES) mesh **is** the surface;
  charges are projected onto mesh vertices with distance weighting.

## 5. The Zwitterionic Metric (Core)

For any surface region *R* with `n⁺` positive and `n⁻` negative charges over
surface area `A`, three components:

- **Balance** (neutrality):
  `B = 1 − |n⁺ − n⁻| / (n⁺ + n⁻)` — 1 when perfectly balanced, 0 when one-sided.
- **Density** (charge-rich, not charge-free):
  `D = (n⁺ + n⁻) / A`, normalized to `D̂ ∈ [0,1]` across the dataset.
- **Mixing** (interspersed, not segregated *within R*): see per-engine below.

**Composite score (product form):**

```
Z(R) = B · D̂ · M
```

The **product** is intentional: any component near zero zeroes the score, which
enforces "balanced **and** dense **and** mixed simultaneously" — exactly the
"neutral but charge-rich and intermixed" definition. (A weighted sum was
rejected because it lets components compensate for one another.)

### Engine A — composite index (prototype)
- Region = Euclidean sphere (radius r ≈ 10 Å, swept) around each surface charge
  center.
- **Mixing** `M_A` = fraction of each charge's short-range neighbors that are of
  **opposite** sign, averaged over R. High when + and − alternate; low when
  like-charges cluster.
- Per-residue Z = aggregate (max/mean) of patch scores covering that residue.
- *Caveat:* Euclidean spheres may bridge across clefts; `M_A` is a coarse
  proxy. Acceptable for the prototype.

### Engine B — bivariate spatial statistics (production)
Treat + and − charges as two **marked point patterns on the SES mesh**, with
distances measured **geodesically** (along the surface, so a patch cannot
"cheat" across a cleft).
- **Mixing** `M_B` derived from the **bivariate pair-correlation / cross-K
  function** `g₊₋(d)`: zwitterionic surfaces show + and − **co-locating at
  short range** (`g₊₋(d) > 1` for small d); segregated surfaces show
  `g₊₋(d) < 1`. This is the formal, literature-grounded "intermixed" measure.
- Computed in **geodesic patches** (fast-marching) so density and balance are
  surface-area-correct.
- Output: a continuous **Z field over mesh vertices**.

**Optional later refinement:** charge **autocorrelation / Moran's I** on the
signed-charge field — true alternation yields *negative* short-range
autocorrelation, an independent confirmation of the mixing signal.

## 6. Patches & Per-Protein Aggregation

**Patch map (from the Z field):**
- Threshold per-location Z (adaptive, e.g., dataset percentile; swept in
  robustness), then **cluster** contiguous high-Z locations into discrete
  patches — connected components on the residue graph (A) or geodesic mesh
  adjacency (B).
- Each patch records: area, mean/max Z, charge composition (n⁺/n⁻), centroid.

**Per-protein score vector** (size/area-normalized so proteins are comparable):
- `Z_frac` — fraction of total surface area that is zwitterionic (coverage).
- `Z_max` — strength of the single best patch.
- `Z_mean` — area-weighted mean Z over the whole surface.
- `n_patches` — count of distinct zwitterionic patches.

This vector is the **Phase-2 seam**: a stability signal attaches as
`corr(stability, Z_*)` with no upstream changes.

## 7. Validation (Phase 1)

Two validation tracks, chosen to establish correctness *without* stability data.

**Positive/negative controls** — a labeled benchmark set:
- *Positives:* halophilic/extremophile proteins (known acidic+basic
  charge-rich surfaces); plus **synthetic positives** (a real surface
  computationally relabeled to perfectly alternating +/−) as an upper-bound
  sanity check.
- *Negatives:* charge-segregated surfaces (strongly dipolar proteins, basic
  DNA-binding patches) and charge-poor/hydrophobic surfaces.
- *Success criterion (set up front):* the metric ranks positives above
  negatives — report AUROC / separation.

**Robustness / sensitivity sweeps** — confirm metric and control ranking are
stable under:
- rSASA threshold, patch radius / geodesic patch size, Z threshold
- pLDDT filtering level
- His fractional weight (0 → 0.1)
- side-chain rotamer noise (small perturbations)
- A-vs-B engine agreement (do the two engines' per-protein scores correlate?)
- *Success criterion (set up front):* control ranking and per-protein score
  *rankings* are stable (Spearman ρ above a pre-set bar) across reasonable
  parameter ranges.

## 8. Tech Stack & Project Layout

**Python**, scientific-bio stack:
- Structures / SASA: `biotite` or `Biopython` + `FreeSASA`
- Mesh (B): `NanoShaper` / `MSMS` → `trimesh` / `pyvista`; geodesics via
  `potpourri3d` / `pygeodesic`
- Spatial stats (B): `pointpats` or a small custom geodesic cross-K function
- Numerics: `numpy` / `scipy` / `pandas`
- Viz: `pyvista` / `matplotlib`
- Config-driven runs for parameter sweeps
- Tests: `pytest` (synthetic surfaces provide exact ground truth for unit tests)

```
prosurf/
  data/            raw structures, control set, cache
  prosurf/
    io/            fetch + parse AF2
    surface/       SASA, charges, mesh (engine A & B)
    metric/        B, D, M, Z  (shared + per-engine)
    patches/       clustering + per-protein aggregation
    validate/      controls, robustness sweeps
    viz/           patch-map rendering
  tests/
  configs/
  docs/superpowers/specs/
```

## 9. Deliverables (Phase 1)

1. Working pipeline with both engines (A prototype, B production) behind a
   common interface.
2. Per-protein Z score vectors + per-protein patch maps for the pilot set.
3. Validation report: control AUROC/separation + robustness sweep results.
4. Rendered patch-map visualizations for representative proteins.

## 10. Phase 2 (out of scope, noted for continuity)

Attach a stability signal (experimental Tm, computational ΔG, pLDDT proxy,
and/or in-cell turnover) at the aggregation seam (§6) and test
`corr(stability, Z_*)` across the pilot and, later, the full human proteome.
