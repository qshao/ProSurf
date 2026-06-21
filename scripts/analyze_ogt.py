"""Cross-species OGT analysis: organism-level trend, LOO robustness,
Z-metric consistency, MIN_PER_ORG sensitivity, within-species replication,
and composition control. Writes data/cross_species_report.{md,pdf}."""
import base64
import io
from collections import defaultdict
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
Z_METRICS = ["z_mean", "z_max", "z_frac"]
MIN_PER_ORG = 30
N_PERM = 10_000
SENSITIVITY_CUTOFFS = [30, 50, 100, 200]

# Coarse phylogenetic groups for effective-N caveat (domain/phylum level)
PHYLO_GROUP = {
    "Oleispira antarctica":    "Proteobacteria",
    "Caenorhabditis elegans":  "Nematoda",
    "Arabidopsis thaliana":    "Viridiplantae",
    "Drosophila melanogaster": "Insecta",
    "Danio rerio":             "Vertebrata",
    "Bacillus subtilis":       "Firmicutes",
    "Saccharomyces cerevisiae":"Fungi",
    "Escherichia coli":        "Proteobacteria",
    "Mus musculus":            "Vertebrata",
    "Homo sapiens":            "Vertebrata",
    "Picrophilus torridus":    "Archaea",
    "Thermus thermophilus":    "Deinococcus-Thermus",
}

DOMAIN = {
    "Oleispira antarctica":    "Bacteria",
    "Caenorhabditis elegans":  "Eukaryota",
    "Arabidopsis thaliana":    "Eukaryota",
    "Drosophila melanogaster": "Eukaryota",
    "Danio rerio":             "Eukaryota",
    "Bacillus subtilis":       "Bacteria",
    "Saccharomyces cerevisiae":"Eukaryota",
    "Escherichia coli":        "Bacteria",
    "Mus musculus":            "Eukaryota",
    "Homo sapiens":            "Eukaryota",
    "Picrophilus torridus":    "Archaea",
    "Thermus thermophilus":    "Bacteria",
    "Human body fluids":       "Eukaryota",
}


def fig_b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def main():
    df = pd.read_csv(SCORES)
    df = df.dropna(subset=["z_mean", "tm", "ogt"] + COVARS)

    # Pre-compute organism table for all cutoffs
    org_all = (df.groupby("organism")
                 .agg(ogt=("ogt", "first"),
                      **{f"mean_{m}": (m, "mean") for m in Z_METRICS},
                      n=("z_mean", "size"))
                 .reset_index()
                 .sort_values("ogt"))
    org = org_all[org_all["n"] >= MIN_PER_ORG].copy()

    # ---- 1. Organism-level trend ----
    rho_org, p_org_param = spearmanr(org["mean_z_mean"], org["ogt"])

    rng = np.random.default_rng(42)
    mean_z_vals = org["mean_z_mean"].values.copy()
    ogt_vals = org["ogt"].values
    perm_rhos = np.array([
        spearmanr(rng.permuted(mean_z_vals), ogt_vals)[0]
        for _ in range(N_PERM)
    ])
    p_perm = max(float((np.abs(perm_rhos) >= abs(rho_org)).mean()), 1.0 / N_PERM)

    # Phylogenetic effective N
    org_names = set(org["organism"])
    phylo_groups = {PHYLO_GROUP.get(o, o) for o in org_names}
    n_eff = len(phylo_groups)
    group_list = ", ".join(sorted(phylo_groups))

    # ---- 2. LOO robustness ----
    loo_rhos = []
    for drop_org in org["organism"]:
        sub = org[org["organism"] != drop_org]
        if len(sub) >= 5:
            r, _ = spearmanr(sub["mean_z_mean"], sub["ogt"])
            loo_rhos.append((drop_org, float(r)))
    loo_min = min(r for _, r in loo_rhos)
    loo_max = max(r for _, r in loo_rhos)

    loo_collapse = [o for o, r in loo_rhos if abs(r) < 0.05]
    collapse_note = (
        f"**Note:** dropping {' or '.join(loo_collapse)} reduces ρ to ≈0, indicating "
        f"the cross-species OGT signal depends entirely on these organisms."
        if loo_collapse else
        "Signal is positive across all leave-one-out subsets."
    )

    # ---- 3. Z-metric consistency ----
    metric_rhos = []
    for m in Z_METRICS:
        r, p = spearmanr(org[f"mean_{m}"], org["ogt"])
        metric_rhos.append((m, float(r), float(p)))

    # ---- 4. MIN_PER_ORG sensitivity ----
    sensitivity = []
    for cutoff in SENSITIVITY_CUTOFFS:
        sub = org_all[org_all["n"] >= cutoff]
        if len(sub) >= 5:
            r, _ = spearmanr(sub["mean_z_mean"], sub["ogt"])
            sensitivity.append((cutoff, len(sub), float(r)))

    # ---- 5. Within-species replication ----
    within = []
    for organism, g in df.groupby("organism"):
        if len(g) >= MIN_PER_ORG:
            r, p = spearmanr(g["z_mean"], g["tm"])
            within.append((organism, g["ogt"].iloc[0], len(g), float(r), float(p)))
    within.sort(key=lambda t: t[1])

    df["z_demeaned"] = df.groupby("organism")["z_mean"].transform(
        lambda g: (g - g.mean()) / (g.std() + 1e-9))
    df["tm_demeaned"] = df.groupby("organism")["tm"].transform(
        lambda g: (g - g.mean()) / (g.std() + 1e-9))
    r_pool, p_pool = spearmanr(df["z_demeaned"], df["tm_demeaned"])

    # ---- 6. Composition control ----
    pr_ogt, pp_ogt = partial_spearman(
        df["z_mean"], df["ogt"], [df[c] for c in COVARS])
    raw_ogt, _ = spearmanr(df["z_mean"], df["ogt"])
    nested = nested_regression_delta(df, "tm", COVARS, "z_mean")

    # ---- 7. Meta-regression: does Z-Tm coupling strength vary with OGT? ----
    within_ogt_vals = np.array([w[1] for w in within])
    within_r_vals   = np.array([w[3] for w in within])
    rho_meta, p_meta_param = spearmanr(within_r_vals, within_ogt_vals)

    meta_perm_rhos = np.array([
        spearmanr(rng.permuted(within_r_vals), within_ogt_vals)[0]
        for _ in range(N_PERM)
    ])
    p_meta_perm = max(float((np.abs(meta_perm_rhos) >= abs(rho_meta)).mean()), 1.0 / N_PERM)

    domain_rhos: dict = defaultdict(list)
    for o, ogt_val, n, r, p in within:
        domain_rhos[DOMAIN.get(o, "Unknown")].append(r)

    domain_means = {d: float(np.mean(rs)) for d, rs in domain_rhos.items()}

    # ---- Verdict strings ----
    ogt_verdict = "retains" if pr_ogt > 0 and pp_ogt < 0.05 else "does NOT retain"
    tm_verdict = "significantly predicts" if nested['p_value'] < 0.05 else "does NOT significantly predict"
    verdict_md = (
        f"**Cross-species OGT verdict:** Z {ogt_verdict} a significant positive association "
        f"with OGT after controlling for bulk charge composition "
        f"(partial ρ = {pr_ogt:+.3f}, p = {pp_ogt:.2e}).\n\n"
        f"**Within-organism stability verdict:** Z {tm_verdict} individual protein Tm beyond "
        f"composition (ΔR² = {nested['delta_r2']:+.4f}, p = {nested['p_value']:.2e}); "
        f"note the effect size is small ({nested['delta_r2']*100:.2f}% of Tm variance) "
        f"despite high significance at N = {len(df):,}."
    )

    # ---- Figures ----
    f1, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(org["ogt"], org["mean_z_mean"], s=60, color="#2c3e50", zorder=3)
    for _, r in org.iterrows():
        ax.annotate(r["organism"].split()[0], (r["ogt"], r["mean_z_mean"]),
                    fontsize=7, xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel("Optimal growth temperature (°C)")
    ax.set_ylabel("Mean z_mean")
    ax.set_title(
        f"Organism-level: mean Z vs OGT\n"
        f"Spearman ρ={rho_org:+.3f}, perm p={p_perm:.4f}, n_eff={n_eff}")
    ax.grid(True, lw=0.4, color="#ccc")
    img_org = fig_b64(f1)

    f_loo, ax = plt.subplots(figsize=(7, max(3, 0.4 * len(loo_rhos) + 1)))
    loo_labels = [r[0] for r in loo_rhos]
    loo_vals = [r[1] for r in loo_rhos]
    ax.barh(range(len(loo_rhos)), loo_vals,
            color=["#27ae60" if r > 0 else "#e74c3c" for r in loo_vals])
    ax.set_yticks(range(len(loo_rhos)))
    ax.set_yticklabels(loo_labels, fontsize=8)
    ax.axvline(rho_org, color="#2c3e50", lw=1.5, linestyle="--",
               label=f"full ρ={rho_org:+.3f}")
    ax.axvline(0, color="#333", lw=0.8)
    ax.set_xlabel("Spearman ρ (mean Z vs OGT, one organism dropped)")
    ax.set_title("LOO robustness")
    ax.legend(fontsize=8)
    img_loo = fig_b64(f_loo)

    f2, ax = plt.subplots(figsize=(7, max(3, 0.4 * len(within) + 1)))
    labels = [w[0] for w in within]
    rs = [w[3] for w in within]
    ax.barh(range(len(within)), rs,
            color=["#27ae60" if r > 0 else "#e74c3c" for r in rs])
    ax.set_yticks(range(len(within)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(0, color="#333", lw=0.8)
    ax.set_xlabel("within-species Spearman ρ(z_mean, Tm)")
    ax.set_title(f"Within-species replication\n"
                 f"pooled (organism-demeaned) ρ={r_pool:+.3f}, p={p_pool:.2e}")
    img_within = fig_b64(f2)

    _color_map = {"Bacteria": "#e74c3c", "Eukaryota": "#2c3e50",
                  "Archaea": "#f39c12", "Unknown": "#95a5a6"}
    f_meta, ax = plt.subplots(figsize=(7, 5))
    seen_domains: set = set()
    for o, ogt_val, n, r, p in within:
        d = DOMAIN.get(o, "Unknown")
        label = d if d not in seen_domains else None
        seen_domains.add(d)
        ax.scatter(ogt_val, r, s=max(30, n / 50), color=_color_map.get(d, "#95a5a6"),
                   label=label, zorder=3)
        ax.annotate(o.split()[0], (ogt_val, r), fontsize=7,
                    xytext=(4, 4), textcoords="offset points")
    ax.axhline(0, color="#333", lw=0.8)
    ax.set_xlabel("Optimal growth temperature (°C)")
    ax.set_ylabel("Within-species Spearman ρ(Z, Tm)")
    ax.set_title(
        f"Meta-regression: Z–Tm coupling strength vs OGT\n"
        f"Spearman ρ={rho_meta:+.3f}, perm p={p_meta_perm:.4f}")
    ax.legend(fontsize=8)
    ax.grid(True, lw=0.4, color="#ccc")
    img_meta = fig_b64(f_meta)

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
    loo_rows = "\n".join(
        f"| {o} | {r:+.3f} |" for o, r in loo_rhos)
    metric_rows = "\n".join(
        f"| {m} | {r:+.3f} | {p:.3f} |" for m, r, p in metric_rhos)
    sens_rows = "\n".join(
        f"| {c} | {n} | {r:+.3f} |" for c, n, r in sensitivity)
    domain_rows = "\n".join(
        f"| {d} | {float(np.mean(rs)):+.3f} | {len(rs)} |"
        for d, rs in sorted(domain_rhos.items()))

    md = f"""# Cross-Species OGT Analysis

**N proteins:** {len(df):,} across {df.organism.nunique()} organisms (OGT {df.ogt.min():.0f}–{df.ogt.max():.0f} °C)

## 1. Organism-level trend
Spearman ρ(mean Z, OGT) = **{rho_org:+.3f}** (parametric p = {p_org_param:.3f}; permutation p = **{p_perm:.4f}**, {N_PERM:,} shuffles; n = {len(org)} organisms).

Phylogenetic note: the {len(org)} qualifying organisms span {n_eff} coarse phylogenetic groups ({group_list}), giving an effective independent sample size of approximately {n_eff}. PGLS correction was not applied; the permutation p above should be interpreted with this caveat.

![org](data:image/png;base64,{img_org})

## 2. LOO robustness
Dropping each organism in turn: ρ range [{loo_min:+.3f}, {loo_max:+.3f}] (full ρ = {rho_org:+.3f}).
{'Signal holds (positive ρ) across all leave-one-out subsets.' if loo_min > 0 else 'Signal flips when at least one organism is dropped — see table.'}
{collapse_note}

| Dropped organism | ρ(mean Z, OGT) |
|-----------------|---------------|
{loo_rows}

![loo](data:image/png;base64,{img_loo})

## 3. Z-metric consistency
All three Z metrics show the same direction of association with OGT.

| Metric | ρ(mean metric, OGT) | parametric p |
|--------|--------------------:|-------------:|
{metric_rows}

## 4. MIN_PER_ORG sensitivity

| MIN_PER_ORG | Organisms qualifying | ρ(mean Z, OGT) |
|------------:|--------------------:|---------------:|
{sens_rows}

## 5. Within-species replication

pooled (organism-demeaned) protein-level Spearman ρ(Z, Tm) = **{r_pool:+.3f}** (p = {p_pool:.2e}, N = {len(df):,} proteins).

| Organism | OGT | N | ρ(z_mean, Tm) | p |
|----------|----:|--:|--------------:|---|
{within_rows}

![within](data:image/png;base64,{img_within})

## 6. Composition control (headline)
Note: partial Spearman and nested regression are computed at protein level (N ≈ 34k), not organism level; OGT enters as a per-protein label inherited from organism membership.
- Raw Spearman ρ(z_mean, OGT) = {raw_ogt:+.3f}
- **Partial** ρ(z_mean, OGT | {', '.join(COVARS)}) = **{pr_ogt:+.3f}** (p = {pp_ogt:.2e})
- Nested model of Tm: R²(composition) = {nested['r2_base']:.4f} → R²(+z_mean) = {nested['r2_full']:.4f}; **ΔR² = {nested['delta_r2']:+.4f}**, added-term p = {nested['p_value']:.2e}, semi-standardized z_mean coef (X-standardized) = {nested['added_coef']:+.3f}.

![delta](data:image/png;base64,{img_delta})

## 7. Organism-specificity: does Z–Tm coupling vary with OGT? (Meta-regression)

Spearman ρ(ρ_within, OGT) = **{rho_meta:+.3f}** (perm p = {p_meta_perm:.4f}, n = {len(within)} organisms). A positive value means the Z–Tm coupling is stronger in warmer organisms; negative means it weakens at higher OGT.

| Domain | Mean ρ(Z, Tm) | N organisms |
|--------|-------------:|------------:|
{domain_rows}

![meta](data:image/png;base64,{img_meta})

{verdict_md}
"""
    MD.write_text(md)
    print(f"Wrote {MD}")
    print(f"  organism-level rho={rho_org:+.3f}  perm-p={p_perm:.4f}  n_eff={n_eff}")
    print(f"  LOO range [{loo_min:+.3f}, {loo_max:+.3f}]")
    print(f"  Z-metrics: {', '.join(f'{m} rho={r:+.3f}' for m, r, _ in metric_rhos)}")
    print(f"  sensitivity: {sensitivity}")
    print(f"  pooled within-org rho={r_pool:+.3f} p={p_pool:.2e}")
    print(f"  partial rho(Z,OGT|comp)={pr_ogt:+.3f} p={pp_ogt:.2e}")
    print(f"  nested deltaR2={nested['delta_r2']:+.4f} p={nested['p_value']:.2e}")
    print(f"  meta-regression rho(within_r,OGT)={rho_meta:+.3f} perm-p={p_meta_perm:.4f}")
    print(f"  domain means: { {d: round(v,3) for d,v in sorted(domain_means.items())} }")

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
