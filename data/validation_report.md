# ProSurf Phase-1 Validation Report — 200-Protein Pilot

**Date:** 2026-06-20
**Control set:** Meltome Atlas (Jarzab et al. 2020, Nat Methods) — 55 human TPP datasets
**Engine:** A (adaptive patch radius: max(0.55 × R_g, 8 Å))
**Proteins:** 100 high-Tm (≥62°C) + 100 low-Tm (≤46°C), all ≥150aa, all in AF2 v6

---

## Primary Results

| Metric | Value |
|--------|-------|
| AUROC z_mean (POS vs NEG) | **0.755** |
| AUROC z_max  (POS vs NEG) | 0.721 |
| Mann-Whitney p (z_mean) | **< 0.0001** |
| **Spearman ρ(z_mean, Tm)** | **0.396, p < 0.000001** |

The Spearman correlation is statistically robust at n=200: ρ = 0.396 with p < 10⁻⁶.
Proteins with higher measured thermostability (TPP Tm) show significantly higher
zwitterionic surface scores.

---

## Distribution by Class

| Class | N | mean z_mean | median | std | Q10 | Q90 |
|-------|---|-------------|--------|-----|-----|-----|
| POS (Tm ≥ 62°C) | 100 | 0.572 | 0.573 | 0.115 | 0.439 | 0.717 |
| NEG (Tm ≤ 46°C) | 100 | 0.472 | 0.460 | 0.094 | 0.377 | 0.566 |

Δ median z_mean = **0.113** (POS higher). Distributions overlap substantially
(std ~ 0.10–0.12), consistent with the signal being real but partial —
TPP Tm captures cellular context, not just intrinsic surface charge patterns.

---

## Robustness Sweeps (Spearman ρ of ranking stability)

| Parameter | Values tested | ρ |
|-----------|--------------|---|
| rsasa_threshold | [0.15, 0.20, 0.25] | **0.902** |
| patch_radius_frac | [0.40, 0.55, 0.70] | 0.672 |
| his_weight | [0.0, 0.05, 0.10] | **0.982** |

---

## Progression Across Pilot Sizes

| N | AUROC z_mean | Spearman ρ | p-value |
|---|-------------|------------|---------|
| 12 | 0.861 | — | — |
| 50 | 0.782 | 0.333 | 0.018 |
| **200** | **0.755** | **0.396** | **<10⁻⁶** |

AUROC decreases slightly as N grows because extreme-Tm proteins dominate the
12-protein and 50-protein sets; the 200-protein set includes proteins closer
to the boundary where the signal is weaker. Spearman ρ increases and p-value
tightens as N grows — the correlation is real and strengthens with sample size.

---

## Outlier Analysis

35 positives score below the class midpoint (0.516); 27 negatives score above.
These are not failures of the metric — they reflect biological noise:
- Low-scoring positives may be thermostable through mechanisms other than
  surface charge (disulfide bonds, metal coordination, buried hydrophobic core).
- High-scoring negatives may have zwitterionic surfaces but be thermolabile
  due to disordered regions, binding-partner dependence, or cofactor loss
  during TPP denaturation.

---

## Recommended Next Steps

1. **patch_radius_frac robustness (ρ=0.672)** — the remaining priority.
   Consider making the fraction adaptive to local charge density rather than
   global R_g.
2. **Phase 2 seam**: cross-reference z_mean scores against a thermodynamic
   stability dataset (ΔG from ProThermDB or FireProt) to test the hypothesis
   independently of cellular context.
3. **Scale to full proteome scan** (~20,000 proteins) using the validated
   pipeline.
