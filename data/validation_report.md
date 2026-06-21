# ProSurf Phase-1 Validation Report

**Date:** 2026-06-20  
**Control set:** Meltome Atlas (Jarzab et al. 2020, Nat Methods) — 55 human TPP datasets  
**Engine:** A (Euclidean spheres, mixing_euclidean)

## Control Set

| UniProt | Gene | Tm (°C) | Class | z_max | z_frac | z_mean | n_patches |
|---------|------|---------|-------|-------|--------|--------|-----------|
| P29401 | TKT | 66.0 | POS | 0.6667 | 0.0462 | 0.6667 | 9 |
| Q96M27 | PRRC1 | 66.8 | POS | 0.8000 | 0.0170 | 0.8000 | 3 |
| Q96GD0 | PDXP | 66.1 | POS | 0.8571 | 0.0414 | 0.6190 | 2 |
| P17174 | GOT1 | 66.8 | POS | 0.8000 | 0.0963 | 0.5333 | 7 |
| Q8NCC3 | PLA2G15 | 67.8 | POS | 1.0000 | 0.0725 | 0.4762 | 6 |
| Q9H6Z4 | RANBP3 | 68.4 | POS | 0.8000 | 0.0442 | 0.4696 | 10 |
| Q92616 | GCN1 | 44.3 | NEG | 0.8571 | 0.0650 | 0.4303 | 44 |
| P32519 | ELF1 | 44.2 | NEG | 0.8000 | 0.0274 | 0.6000 | 7 |
| Q6NXE6 | ARMC6 | 43.2 | NEG | 1.0000 | 0.1283 | 0.4216 | 14 |
| Q8IX90 | SKA3 | 44.0 | NEG | 1.0000 | 0.0681 | 0.4048 | 11 |
| Q15742 | NAB2 | 43.2 | NEG | 0.6667 | 0.0334 | 0.6667 | 6 |
| P04264 | KRT1 | 43.2 | NEG | 0.8000 | 0.0481 | 0.5161 | 17 |

## AUROC (positives = high Tm ≥ 66°C, negatives = low Tm ≤ 44.3°C)

| Metric | AUROC | Interpretation |
|--------|-------|----------------|
| z_max  | 0.417  | Best-patch score — inversely ordered (negatives outscore positives) |
| z_frac | 0.444  | Zwitterionic surface fraction — slightly inverted |
| z_mean | 0.736  | Mean patch score — positives slightly higher |

## Robustness Sweeps (Spearman ρ of protein ranking across parameter values)

| Parameter | Values tested | Spearman ρ |
|-----------|--------------|------------|
| rsasa_threshold | [0.15, 0.20, 0.25] | 0.490 |
| patch_radius    | [8.0, 10.0, 12.0]  | 0.041 |
| his_weight      | [0.0, 0.05, 0.10]  | 0.947 |

## Engine A vs B Agreement

Spearman ρ (z_max): **1.000** (perfect agreement — cross_colocation gives identical z_max on this dataset)

## Key Findings

1. **z_mean AUROC = 0.736** — the strongest signal; positives trend toward higher mean patch Z scores.
2. **z_max and z_frac are inverted** — larger low-Tm proteins (GCN1 2671aa, ARMC6 501aa, SKA3 412aa) accumulate many patches, driving high z_max and z_frac despite low thermostability.
3. **patch_radius robustness is poor (ρ=0.041)** — protein rankings change substantially between 8–12 Å. This is the most actionable metric design issue.
4. **his_weight robustness is excellent (ρ=0.947)** — excluding His does not affect rankings.
5. **Engine A = Engine B exactly** on this dataset — cross_colocation produces the same z_max as mixing_euclidean (the D̂ normalization dominates; mixing component provides no additional discrimination here).

## Recommended Next Steps

- **Size-normalize the metric**: divide n_patches by surface area or protein length to avoid size confound.
- **Fix patch_radius sensitivity**: the ρ=0.041 result means the ranking is essentially random across 8–12 Å. Likely cause: patch_radius controls cluster adjacency but many patches change identity when radius changes. Consider making it adaptive (% of protein radius) rather than absolute.
- **z_mean as primary metric**: shows the best discriminative signal in this pilot.
- **Expand to 50–100 proteins** with Tm data to get stable AUROC estimates.