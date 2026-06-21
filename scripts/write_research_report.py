#!/usr/bin/env python3
"""
Generate a comprehensive scientific research report for the ProSurf cross-species OGT analysis.

Outputs:
  data/ogt_research_report.md   — full markdown report with embedded figures
  data/ogt_research_report.pdf  — PDF render via weasyprint + markdown
"""

import base64
import io
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse, FancyBboxPatch
import numpy as np

# ── paths ────────────────────────────────────────────────────────────────────
REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
SOURCE_REPORT = DATA / "cross_species_report.md"
OUT_MD = DATA / "ogt_research_report.md"
OUT_PDF = DATA / "ogt_research_report.pdf"


# ── helper: encode figure to base64 PNG ──────────────────────────────────────
def fig_b64(fig: plt.Figure) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("ascii")


# ── Figure 1: conceptual / motivation figure ─────────────────────────────────
def make_concept_figure() -> str:
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(14, 6))

    # ── Panel A — Zwitterionic Surface Score Z ────────────────────────────
    ax_a.set_xlim(0, 10)
    ax_a.set_ylim(0, 10)
    ax_a.set_aspect("equal")
    ax_a.axis("off")
    ax_a.set_title("(A) Zwitterionic Surface Score", fontsize=12, fontweight="bold", pad=10)

    # protein ellipse
    ellipse = Ellipse((5, 5), 7, 5, facecolor="#d5d8dc", edgecolor="#7f8c8d",
                      linewidth=2, zorder=1)
    ax_a.add_patch(ellipse)

    # 12 charge labels alternating +/−
    angles = np.linspace(0, 2 * np.pi, 12, endpoint=False)
    signs = ["+", "−"] * 6
    colors = ["#e74c3c", "#2980b9"] * 6
    for theta, sign, col in zip(angles, signs, colors):
        x = 3.5 * np.cos(theta) + 5
        y = 2.5 * np.sin(theta) + 5
        ax_a.text(x, y, sign, ha="center", va="center", fontsize=14,
                  fontweight="bold", color=col, zorder=2)

    # formula and legend below ellipse
    ax_a.text(5, 1.8, r"$Z = B \cdot \hat{D} \cdot M$",
              ha="center", va="center", fontsize=14, fontweight="bold")
    ax_a.text(5, 1.1, "Balance × Density × Mixing",
              ha="center", va="center", fontsize=9, color="#7f8c8d")

    # ── Panel B — Research Workflow ───────────────────────────────────────
    ax_b.set_xlim(0, 10)
    ax_b.set_ylim(0, 10)
    ax_b.axis("off")
    ax_b.set_title("(B) Research Workflow", fontsize=12, fontweight="bold", pad=10)

    boxes = [
        (8.0, "#3498db", "Meltome Atlas\n13 organisms · OGT 15–70°C\nPer-protein Tm (thermal proteomics)"),
        (6.2, "#27ae60", "AlphaFold2 Structures\n34,231 proteins via EBI API\nLength filter ≥ 150 aa"),
        (4.4, "#e67e22", r"ProSurf Scoring" + "\nZ = B · D̂ · M per protein\n+ Composition control features"),
        (2.6, "#9b59b6", "Statistical Analysis\nOrganism-level Spearman · LOO\nPartial correlation · Meta-regression"),
    ]

    box_h = 1.0
    box_w = 5.0
    box_x = 5.0  # centre x

    box_tops = []
    box_bottoms = []
    for y_centre, colour, text in boxes:
        left = box_x - box_w / 2
        bottom = y_centre - box_h / 2
        patch = FancyBboxPatch(
            (left, bottom), box_w, box_h,
            boxstyle="round,pad=0.1",
            facecolor=colour, edgecolor="none", zorder=2
        )
        ax_b.add_patch(patch)
        ax_b.text(box_x, y_centre, text,
                  ha="center", va="center", fontsize=9,
                  fontweight="bold", color="white", zorder=3,
                  multialignment="center")
        box_tops.append(y_centre + box_h / 2)
        box_bottoms.append(y_centre - box_h / 2)

    # arrows between boxes
    for i in range(len(box_tops) - 1):
        y_start = box_bottoms[i]
        y_end = box_tops[i + 1]
        ax_b.annotate(
            "", xy=(box_x, y_end), xytext=(box_x, y_start),
            arrowprops=dict(arrowstyle="->", color="#2c3e50", lw=1.5),
            zorder=1
        )

    plt.tight_layout()
    encoded = fig_b64(fig)
    plt.close(fig)
    return encoded


# ── Extract existing figures from cross_species_report.md ────────────────────
def extract_figures(report_path: Path) -> dict:
    text = report_path.read_text(encoding="utf-8")
    pattern = r'!\[([^\]]*)\]\(data:image/png;base64,([A-Za-z0-9+/=]+)\)'
    figures = {}
    for m in re.finditer(pattern, text):
        tag = m.group(1)
        b64 = m.group(2)
        figures[tag] = b64
    return figures


# ── Report template ───────────────────────────────────────────────────────────
REPORT_TEMPLATE = r"""# Zwitterionic Surface Charge Mixing as a Predictor of Protein Thermal Stability: A Cross-Species Analysis

---

## Abstract

Protein thermal stability is essential for enzymatic function across the temperature range of life. The ProSurf zwitterionic surface score Z — quantifying the balance, density, and spatial mixing of surface charged residues — has been proposed as a structural correlate of thermostability. Here we test whether Z rises with optimal growth temperature (OGT) across the tree of life using 34,231 AlphaFold2-predicted protein structures from 13 organisms spanning OGT 15–70°C, with melting temperatures from the Meltome Atlas [1]. Organism-level Z does not significantly correlate with OGT (Spearman ρ = +0.21, permutation p = 0.50, n = 13 organisms), and this trend vanishes after controlling for bulk charge composition (partial ρ = +0.006, p = 0.30). However, within individual organisms, proteins with higher melting temperature (Tm) consistently show higher Z (pooled within-organism Spearman ρ = +0.079, p = 1.26×10⁻⁴⁸, N = 34,231), and Z significantly improves Tm prediction beyond composition features (ΔR² = +0.0038, p = 3.67×10⁻³⁰). This within-organism effect is strongest in bacteria (mean ρ = +0.096) and absent in the single archaeal representative (ρ = +0.031). Zwitterionic surface organisation is therefore a genuine but modest predictor of individual protein stability, rather than a driver of cross-species thermal adaptation.

---

## 1. Introduction

Life on Earth spans an extraordinary temperature range, from psychrophiles active near 0°C to hyperthermophiles growing above 100°C [4]. Proteins must maintain structural integrity and catalytic function across these temperatures, imposing strong selective pressure on thermostability. Known molecular strategies include increased hydrophobic core packing, reduced loop length, introduction of disulfide bridges, and electrostatic stabilisation through salt bridge networks [4, 5]. How surface charge *organisation* — distinct from overall charge content — contributes to stability has received less systematic attention at a proteome scale.

Zwitterionic materials, characterised by intimately paired positive and negative charges, are known to exhibit unusual interfacial properties including enhanced hydration and resistance to non-specific interactions [6]. By analogy, protein surfaces with balanced, densely mixed opposite charges may form local electrostatic networks that resist thermal unfolding. We term this property the **zwitterionic surface score Z**, computed as the product of three geometric components: Balance (B, equality of positive/negative surface charge counts), normalised Density (D̂, charged residue density on the solvent-exposed surface), and Mixing (M, spatial interleaving of opposite charges at the residue level). Together, Z = B · D̂ · M captures the degree to which a protein surface presents an evenly distributed, spatially mixed charge pattern (Figure 1A).

The Meltome Atlas [1] provides an exceptional resource to test the hypothesis that Z increases with OGT: thermal proteome profiling data for 13 organisms spanning OGT 15–70°C, yielding per-protein Tm values for tens of thousands of proteins. Structural predictions from AlphaFold2 [2, 3] enable proteome-scale Z computation without dependence on experimentally solved structures.

We test two nested hypotheses (Figure 1B): (H1) organisms adapted to higher OGT have proteomes with higher mean Z; and (H2) within individual organisms, thermally stable proteins have higher Z than labile ones. We additionally ask whether any Z-OGT or Z-Tm relationship is independent of bulk charge composition — a critical control because thermophilic proteomes tend to be enriched in charged residues overall [5].

![Figure 1: Motivation and research design](data:image/png;base64,INSERT_CONCEPT)

**Figure 1.** (A) The ProSurf zwitterionic surface score Z = B · D̂ · M captures three geometric aspects of surface charge organisation: Balance (equality of positive and negative charges), normalised Density (charged residue surface density), and Mixing (spatial interleaving of opposite charges). (B) Research workflow: Meltome Atlas thermal proteomics data provide per-protein Tm values and organism OGTs; AlphaFold2 structures (EBI API) provide atomic coordinates; ProSurf computes Z and composition features; statistical analysis tests H1 (cross-species OGT trend) and H2 (within-organism Tm association).

---

## 2. Methods

### 2.1 Meltome Atlas data

Melting temperatures and organism metadata were obtained from the Meltome Atlas (Jarzab et al., 2020) [1], which reports per-protein Tm values from thermal proteome profiling experiments. The per-organism optimal growth temperatures were taken from Supplementary Table S1. For proteins measured across multiple replicate datasets from the same organism, the median Tm was used. We included only canonical UniProt accessions matching the pattern `[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]` or `[OPQ][0-9][A-Z0-9]{3}[0-9]`.

### 2.2 Protein structures

Atomic coordinates were obtained from the AlphaFold Protein Structure Database [3] via the EBI REST API, requesting AlphaFold2 version 6 models. Proteins shorter than 150 residues were excluded to ensure reliable surface area estimation. Human proteins with pre-existing scores from a full-proteome run were seeded from that dataset; composition features were computed from locally cached PDB files where available.

### 2.3 ProSurf Z score

Surface residues were identified by relative solvent-accessible surface area (rSASA threshold = 0.05, computed using the Shrake–Rupley algorithm [7] as implemented in Biotite [8]). Charged residue assignments: K and R receive charge +1; D and E receive charge −1; H receives +0.1 (his_weight). The three Z components are:

- **Balance** B = min(n₊, n₋) / max(n₊, n₋), where n₊ and n₋ are surface positive and negative charge totals
- **Density** D̂ = (n₊ + n₋) / SASA, normalised to [0, 1] across the dataset
- **Mixing** M = mean, over all surface charged residues, of the fraction of their charged neighbours carrying the opposite sign

Z = B · D̂ · M. Three Z metrics were computed: z_mean (mean patch Z), z_max (maximum patch Z), z_frac (fraction of patches with Z > 0).

### 2.4 Composition control features

Four bulk charge features were computed as controls:
- *charged_frac*: fraction of all residues that are D, E, K, or R
- *surface_charged_frac*: fraction of surface residues that are charged
- *net_charge_per_res*: |n₊ − n₋| / n_residues (surface charges, unnormalised imbalance)
- *surface_charge_density*: (n₊ + n₋) / SASA (Å⁻²)

These features capture the *amount* and *imbalance* of surface charge but not its spatial mixing.

### 2.5 Statistical analysis

**Organism-level:** Spearman correlation between per-organism mean Z and OGT. A permutation p-value (10,000 shuffles of organism-level mean Z; floor at 1/10,000) replaced the parametric p, which is unreliable at N ≤ 15 organisms.

**Leave-one-out (LOO):** Each of the 13 organisms was dropped in turn; the Spearman ρ was recomputed on the remaining subset (minimum 5 organisms required).

**Phylogenetic effective N:** Organisms were mapped to coarse phylogenetic groups (domain/phylum level), and the number of distinct groups was taken as an approximate effective independent sample size.

**Partial Spearman:** The Z-OGT association was recomputed after regressing out composition covariates (COVARS = [charged_frac, net_charge_per_res, surface_charge_density]) from both Z and OGT by OLS on rank-transformed variables.

**Nested OLS regression:** Two models of protein Tm were compared: (base) composition covariates only; (full) composition + z_mean. Predictors were z-standardised. The incremental ΔR² was tested with an F-test.

**Pooled within-organism Spearman:** Z and Tm were separately demeaned and z-scored within each organism; the resulting residuals were pooled and Spearman-correlated.

**Meta-regression:** The within-organism Spearman ρ(Z, Tm) for each organism was correlated with that organism's OGT (Spearman, with 10,000-shuffle permutation p). Organisms were grouped into Bacteria, Eukaryota, and Archaea; mean ρ was computed per group.

All analyses required ≥ 30 scored proteins per organism (MIN_PER_ORG = 30). Analyses were performed in Python 3.11 using SciPy [9], NumPy [10], pandas [11], and matplotlib [12]. Structures were parsed with Biotite [8].

---

## 3. Results

### 3.1 Dataset overview

After quality filtering (≥ 150 residues, ≥ 30 proteins per organism), 34,231 proteins across 13 organisms were retained, spanning OGT 15°C (*Oleispira antarctica*) to 70°C (*Thermus thermophilus*). The dataset includes psychrophiles, mesophiles, and two thermophiles (*Picrophilus torridus*, OGT 60°C; *Thermus thermophilus*, OGT 70°C). Organism sizes ranged from 100 proteins (Human body fluids dataset) to 10,695 (*Homo sapiens*).

### 3.2 Organism-level Z–OGT trend

The Spearman correlation between per-organism mean Z (z_mean) and OGT was ρ = +0.21 (permutation p = 0.50, n = 13 organisms; Figure 2). All three Z metrics showed positive but non-significant associations (z_mean ρ = +0.21, z_max ρ = +0.16, z_frac ρ = +0.15; Table 1). The result was stable across minimum protein-count thresholds of 30–200 (Table 2). The 13 organisms span 10 coarse phylogenetic groups, giving an effective independent sample size of approximately 10; phylogenetic non-independence may further reduce effective power.

![Figure 2: Organism-level Z vs OGT](data:image/png;base64,INSERT_ORG)

**Figure 2.** Organism-level mean Z (z_mean) vs optimal growth temperature (OGT). Each point represents one organism (N = 13). Spearman ρ = +0.21, permutation p = 0.50 (10,000 shuffles). The trend is in the expected direction but not statistically significant.

**Table 1.** Z-metric consistency: Spearman ρ between per-organism mean metric and OGT.

| Metric | ρ(mean metric, OGT) | Parametric p |
|--------|--------------------:|-------------:|
| z_mean | +0.212 | 0.486 |
| z_max  | +0.162 | 0.597 |
| z_frac | +0.154 | 0.616 |

**Table 2.** Sensitivity to minimum protein-count threshold (MIN_PER_ORG).

| MIN_PER_ORG | Organisms qualifying | ρ(mean Z, OGT) |
|------------:|--------------------:|---------------:|
| 30  | 13 | +0.212 |
| 50  | 13 | +0.212 |
| 100 | 13 | +0.212 |
| 200 | 12 | +0.297 |

### 3.3 LOO robustness

Leave-one-organism-out analysis revealed that the positive trend is entirely driven by the two thermophilic organisms: dropping either *Picrophilus torridus* or *Thermus thermophilus* reduces the correlation to exactly ρ = 0.000 (Figure 3; Table 3). With all mesophilic/psychrophilic organisms alone (OGT 15–37°C), no Z–OGT trend is detectable.

![Figure 3: LOO robustness](data:image/png;base64,INSERT_LOO)

**Figure 3.** Leave-one-organism-out Spearman ρ(mean Z, OGT). Horizontal dashed line shows the full-dataset ρ = +0.212. Dropping either thermophile (*Picrophilus* or *Thermus*) collapses the correlation to zero.

**Table 3.** Leave-one-organism-out results.

| Dropped organism | OGT (°C) | ρ(mean Z, OGT) |
|-----------------|----------:|---------------:|
| *Oleispira antarctica* | 15 | +0.286 |
| *Caenorhabditis elegans* | 20 | +0.279 |
| *Arabidopsis thaliana* | 25 | +0.279 |
| *Danio rerio* | 28 | +0.257 |
| *Drosophila melanogaster* | 28 | +0.193 |
| *Bacillus subtilis* | 30 | +0.175 |
| *Saccharomyces cerevisiae* | 30 | +0.260 |
| *Escherichia coli* | 37 | +0.081 |
| *Homo sapiens* | 37 | +0.329 |
| *Mus musculus* | 37 | +0.297 |
| *Picrophilus torridus* | 60 | +0.000 |
| *Thermus thermophilus* | 70 | +0.000 |

### 3.4 Within-organism Z–Tm association

Within individual organisms, proteins with higher Tm showed consistently higher Z (Figure 4). The pooled within-organism Spearman ρ was +0.079 (p = 1.26×10⁻⁴⁸, N = 34,231 proteins), statistically robust despite the modest effect size.

Per-organism effects were positive and significant in 8 of 13 organisms (Table 4). The strongest associations were in bacteria: *Bacillus subtilis* (ρ = +0.149), *Escherichia coli* (ρ = +0.109), *Thermus thermophilus* (ρ = +0.087). Vertebrate proteomes showed moderate effects (*Homo sapiens* ρ = +0.086, *Mus musculus* ρ = +0.081). *Drosophila melanogaster* (ρ = −0.026) and *Picrophilus torridus* (ρ = +0.031) showed no significant relationship.

![Figure 4: Within-species replication](data:image/png;base64,INSERT_WITHIN)

**Figure 4.** Within-species Spearman ρ(Z, Tm) for each organism, sorted by OGT. Bar colour: green = positive, red = negative. Title reports pooled within-organism ρ across all 34,231 proteins.

**Table 4.** Within-organism Spearman ρ(z_mean, Tm).

| Organism | OGT (°C) | N | ρ(z_mean, Tm) | p |
|----------|----------:|--:|--------------:|---|
| *Oleispira antarctica* | 15 | 1,261 | +0.038 | 0.173 |
| *Caenorhabditis elegans* | 20 | 3,032 | +0.041 | 0.023 |
| *Arabidopsis thaliana* | 25 | 2,380 | +0.098 | 1.69×10⁻⁶ |
| *Danio rerio* | 28 | 1,786 | +0.078 | 1.04×10⁻³ |
| *Drosophila melanogaster* | 28 | 1,560 | −0.026 | 0.310 |
| *Bacillus subtilis* | 30 | 1,354 | +0.149 | 3.88×10⁻⁸ |
| *Saccharomyces cerevisiae* | 30 | 1,757 | +0.086 | 2.87×10⁻⁴ |
| *Escherichia coli* | 37 | 1,987 | +0.109 | 1.14×10⁻⁶ |
| *Homo sapiens* | 37 | 10,695 | +0.086 | 5.45×10⁻¹⁹ |
| *Mus musculus* | 37 | 6,434 | +0.081 | 6.91×10⁻¹¹ |
| *Picrophilus torridus* | 60 | 767 | +0.031 | 0.385 |
| *Thermus thermophilus* | 70 | 1,118 | +0.087 | 3.46×10⁻³ |

### 3.5 Composition control

At protein level, the raw Spearman ρ(z_mean, OGT) = −0.041 (where OGT is a per-protein constant inherited from organism). After controlling for bulk charge composition (charged_frac, net_charge_per_res, surface_charge_density), the partial Spearman ρ = +0.006 (p = 0.30), confirming that Z carries no independent OGT signal beyond bulk charge content.

In the nested OLS model of protein Tm, however, Z adds significantly beyond the composition-only baseline: R²(composition) = 0.0086 → R²(+Z) = 0.0124; ΔR² = +0.0038 (p = 3.67×10⁻³⁰; Figure 5). The semi-standardised Z coefficient is +0.688 °C/σ_Z. The effect size is small — Z accounts for 0.38% of Tm variance — but the signal is independently confirmed by the pooled within-organism Spearman.

![Figure 5: Nested model composition control](data:image/png;base64,INSERT_DELTA)

**Figure 5.** Nested OLS model of protein Tm. Left bar: R² for composition features alone (charged_frac, net_charge_per_res, surface_charge_density). Right bar: R² after adding z_mean. ΔR² = +0.0038, p = 3.67×10⁻³⁰.

### 3.6 Organism-specificity: meta-regression

The meta-regression of within-organism ρ(Z, Tm) against OGT yielded ρ = +0.117 (permutation p = 0.71; Figure 6), indicating no significant trend in Z–Tm coupling strength across the OGT range. However, stratifying by domain of life revealed a consistent gradient: bacteria showed the strongest mean Z–Tm coupling (mean ρ = +0.096, n = 4 organisms), followed by eukaryotes (mean ρ = +0.044, n = 8), with archaea showing the weakest association (mean ρ = +0.031, n = 1).

![Figure 6: Meta-regression](data:image/png;base64,INSERT_META)

**Figure 6.** Meta-regression: within-organism Spearman ρ(Z, Tm) vs OGT for each organism, coloured by domain of life. Point size is proportional to number of proteins scored. Spearman ρ(ρ_within, OGT) = +0.117, permutation p = 0.71.

---

## 4. Discussion

### 4.1 Cross-species OGT trend

The primary hypothesis — that zwitterionic surface character rises systematically with OGT — is not supported at the organism level after composition control. The observed ρ = +0.21 is non-significant even before controlling for composition, and drops to ρ = +0.006 thereafter. Furthermore, the LOO analysis demonstrates that the entire residual positive trend depends on the two most extreme thermophiles in the dataset (*Picrophilus torridus* and *Thermus thermophilus*): removing either reduces the correlation to zero. With only two thermophilic data points and ten effective independent phylogenetic groups, the cross-species test is severely underpowered.

These findings are consistent with the known biochemistry of thermophilic adaptation. Thermophilic proteomes have more charged residues overall [5], and these charges are deployed primarily in salt bridge networks that provide enthalpic stabilisation. Our composition covariates (charged_frac, surface_charge_density) capture this bulk increase, leaving Z — which specifically measures *mixing* and *balance* of surface charges — with no residual OGT signal. In other words, thermophiles appear to increase the *amount* of surface charge rather than the *spatial organisation* of those charges.

### 4.2 Within-organism Z–Tm relationship

In contrast, the within-organism analysis reveals a genuine, reproducible relationship between Z and thermal stability at the protein level (pooled ρ = +0.079, p = 1.26×10⁻⁴⁸). This finding survives composition control (ΔR² = +0.0038, p = 3.67×10⁻³⁰), indicating that the Z–Tm association at the protein level is not simply a proxy for charge-rich proteins being more stable. Proteins within a given proteome that have evolved higher stability tend to have better-mixed surface charge patterns, independent of how many charged residues they carry.

The effect size (0.38% of Tm variance explained) is small, which is unsurprising: protein thermal stability is a complex phenotype determined by many factors including hydrophobic packing, hydrogen bonding, loop entropy, and cofactor binding [4, 5, 7]. Z captures only one aspect of surface organisation. The biological relevance of the 0.38% explained variance would need to be assessed in the context of the full model of Tm determinants.

### 4.3 Domain-of-life specificity

The observation that bacteria show consistently stronger Z–Tm coupling than eukaryotes (mean ρ bacteria = +0.096 vs eukaryotes = +0.044) is intriguing. Several non-exclusive explanations are possible. First, bacterial proteins are on average smaller and structurally simpler than eukaryotic proteins, potentially making surface electrostatics a more dominant determinant of stability relative to other factors. Second, eukaryotic proteins are more commonly stabilised by post-translational modifications (glycosylation, disulfide bonds, phosphorylation) and by chaperone-assisted folding pathways, which may reduce the requirement for intrinsic surface charge organisation. Third, the eukaryotic organisms in our dataset include organisms as divergent as *Arabidopsis*, *Drosophila*, and *Homo sapiens*, with substantially different proteome compositions and cellular environments.

The archaea result (ρ = +0.031, not significant) is based on a single organism (*Picrophilus torridus*) and should not be over-interpreted. *Picrophilus* is an extreme acidothermophile (OGT 60°C, pH optima ~0.7) whose proteins may employ unusual stabilisation mechanisms — including extensive protonation of surface residues — that decouple Z from thermal stability.

### 4.4 Limitations

Several limitations constrain the conclusions. First, the cross-species test is underpowered (N = 13 organisms, n_eff ≈ 10 phylogenetic groups), and phylogenetic correlations among closely related organisms (e.g., three vertebrate representatives) are not formally corrected. PGLS correction [13] would require a calibrated species tree and was not performed. Second, AlphaFold2 structures represent predicted conformations that may not fully capture in-solution dynamics or co-factor-induced remodelling; surface charge environments could differ from the crystallographic or solution state. Third, Tm values from thermal proteome profiling are operationally defined as the temperature at which a protein precipitates from lysate rather than the thermodynamic melting temperature, which may introduce measurement noise. Fourth, the Z score itself is a geometric simplification of a complex electrostatic surface, and alternative charge-mixing metrics might yield different effect sizes.

---

## 5. Conclusion

We performed a systematic cross-species analysis of the ProSurf zwitterionic surface score Z across 34,231 proteins from 13 organisms spanning OGT 15–70°C. The cross-species OGT trend (ρ = +0.21) is not statistically significant and vanishes after controlling for bulk charge composition, with the entire apparent trend attributable to two thermophilic organisms. In contrast, Z is a consistent and composition-independent predictor of individual protein thermal stability within organisms (ΔR² = +0.0038, p = 3.67×10⁻³⁰), most prominently in bacteria. These results indicate that zwitterionic surface charge mixing shapes intra-proteome variation in protein stability but does not account for cross-species thermal adaptation at the OGT level.

---

## References

[1] Jarzab, A., Kurzawa, N., Hopf, T., Moerch, M., Zecha, J., Meng, N., ... & Kuster, B. (2020). Meltome atlas — thermal proteome stability across the tree of life. *Nature Methods*, 17, 495–503.

[2] Jumper, J., Evans, R., Pritzel, A., Green, T., Figurnov, M., Ronneberger, O., ... & Hassabis, D. (2021). Highly accurate protein structure prediction with AlphaFold. *Nature*, 596, 583–589.

[3] Varadi, M., Anyango, S., Deshpande, M., Nair, S., Natassia, C., Yordanova, G., ... & Velankar, S. (2022). AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. *Nucleic Acids Research*, 50, D439–D444.

[4] Sterner, R., & Liebl, W. (2001). Thermophilic adaptation of proteins. *Critical Reviews in Biochemistry and Molecular Biology*, 36, 39–106.

[5] Vieille, C., & Zeikus, G. J. (2001). Hyperthermophilic enzymes: sources, uses, and molecular mechanisms for thermostability. *Microbiology and Molecular Biology Reviews*, 65, 1–43.

[6] Chen, S., Li, L., Zhao, C., & Zheng, J. (2010). Surface hydration: Principles and applications toward low-fouling/nonfouling biomaterials. *Polymer*, 51, 5283–5293.

[7] Shrake, A., & Rupley, J. A. (1973). Environment and exposure to solvent of protein atoms: lysozyme and insulin. *Journal of Molecular Biology*, 79, 351–371.

[8] Kunzmann, P., & Hamacher, K. (2018). Biotite: a unifying open source computational biology framework in Python. *BMC Bioinformatics*, 19, 346.

[9] Virtanen, P., Gommers, R., Oliphant, T. E., Haberland, M., Reddy, T., Cournapeau, D., ... & SciPy 1.0 Contributors (2020). SciPy 1.0: Fundamental algorithms for scientific computing in Python. *Nature Methods*, 17, 261–272.

[10] Harris, C. R., Millman, K. J., van der Walt, S. J., Gommers, R., Virtanen, P., Cournapeau, D., ... & Oliphant, T. E. (2020). Array programming with NumPy. *Nature*, 585, 357–362.

[11] McKinney, W. (2010). Data structures for statistical computing in Python. In *Proceedings of the 9th Python in Science Conference* (Vol. 445, pp. 51–56).

[12] Hunter, J. D. (2007). Matplotlib: A 2D graphics environment. *Computing in Science & Engineering*, 9, 90–95.

[13] Grafen, A. (1989). The phylogenetic regression. *Philosophical Transactions of the Royal Society of London B*, 326, 119–157.
"""


# ── PDF rendering ─────────────────────────────────────────────────────────────
def render_pdf(md_text: str, out_path: Path) -> bool:
    """Render markdown to PDF via weasyprint. Returns True on success."""
    try:
        import markdown as md_lib
        import weasyprint

        html_body = md_lib.markdown(
            md_text,
            extensions=["tables", "fenced_code", "toc"],
        )

        # Minimal CSS for a readable scientific report
        css = """
        @page { size: A4; margin: 2.5cm; }
        body {
            font-family: Georgia, serif;
            font-size: 11pt;
            line-height: 1.6;
            color: #111;
            max-width: 100%;
        }
        h1 { font-size: 16pt; text-align: center; margin-bottom: 0.5em; }
        h2 { font-size: 13pt; border-bottom: 1px solid #ccc; padding-bottom: 3px; margin-top: 1.5em; }
        h3 { font-size: 11pt; margin-top: 1.2em; }
        p  { text-align: justify; margin: 0.5em 0; }
        table {
            border-collapse: collapse;
            width: 100%;
            font-size: 9pt;
            margin: 1em 0;
        }
        th, td { border: 1px solid #aaa; padding: 4px 8px; }
        th { background: #eee; font-weight: bold; }
        img { max-width: 100%; height: auto; display: block; margin: 1em auto; }
        code { font-family: monospace; background: #f4f4f4; padding: 1px 3px; }
        hr { border: none; border-top: 1px solid #ccc; margin: 1.5em 0; }
        """

        full_html = f"""<!DOCTYPE html>
<html lang="en">
<head><meta charset="utf-8"><style>{css}</style></head>
<body>{html_body}</body>
</html>"""

        doc = weasyprint.HTML(string=full_html)
        doc.write_pdf(str(out_path))
        return True
    except Exception as exc:
        print(f"PDF render failed: {exc}", file=sys.stderr)
        return False


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    # 1. Generate conceptual figure
    print("Generating conceptual figure ...", flush=True)
    img_concept = make_concept_figure()

    # 2. Extract existing analysis figures
    print(f"Extracting figures from {SOURCE_REPORT} ...", flush=True)
    figures = extract_figures(SOURCE_REPORT)
    missing = [t for t in ("org", "loo", "within", "delta", "meta") if t not in figures]
    if missing:
        print(f"WARNING: missing figure tags: {missing}", file=sys.stderr)

    # 3. Assemble report
    print("Assembling report ...", flush=True)
    report = REPORT_TEMPLATE
    report = report.replace("INSERT_CONCEPT", img_concept)
    for tag in ("org", "loo", "within", "delta", "meta"):
        placeholder = f"INSERT_{tag.upper()}"
        if tag in figures:
            report = report.replace(placeholder, figures[tag])
        else:
            report = report.replace(placeholder, "")

    # 4. Write markdown
    OUT_MD.write_text(report, encoding="utf-8")
    word_count = len(report.split())
    print(f"Wrote {OUT_MD}")
    print(f"Word count: {word_count:,}")

    # 5. Render PDF
    print("Rendering PDF ...", flush=True)
    pdf_ok = render_pdf(report, OUT_PDF)
    if pdf_ok:
        print(f"Wrote {OUT_PDF}")
    else:
        print("PDF not generated (see stderr for details).")

    # 6. Verify sections
    required_sections = [
        "## Abstract",
        "## 1. Introduction",
        "## 2. Methods",
        "## 3. Results",
        "## 4. Discussion",
        "## 5. Conclusion",
        "## References",
    ]
    for sec in required_sections:
        if sec not in report:
            print(f"WARNING: section missing: {sec}", file=sys.stderr)
        else:
            print(f"  OK  {sec}")


if __name__ == "__main__":
    main()
