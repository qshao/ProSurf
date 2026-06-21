#!/usr/bin/env python3
"""
Generate a summary slide deck PDF for the ProSurf cross-species OGT analysis.

Output: report/ogt_slide_deck.pdf
"""

import base64
import io
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import FancyBboxPatch, Ellipse
import numpy as np

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
REPORT_DIR = REPO / "report"
REPORT_DIR.mkdir(exist_ok=True)
SOURCE_REPORT = DATA / "cross_species_report.md"
OUT_PDF = REPORT_DIR / "ogt_slide_deck.pdf"

W, H = 13.33, 7.5  # 16:9 inches

# ── Palette ───────────────────────────────────────────────────────────────────
NAVY   = "#1a3a5c"
BLUE   = "#2980b9"
GREEN  = "#27ae60"
ORANGE = "#e67e22"
RED    = "#e74c3c"
DARK   = "#2c3e50"
LGREY  = "#f0f4f8"
WHITE  = "white"


# ── Figure extraction ─────────────────────────────────────────────────────────
def extract_figures(path: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    figs = {}
    for m in re.finditer(
        r'!\[([^\]]*)\]\(data:image/png;base64,([A-Za-z0-9+/=\n]+)\)', text
    ):
        tag = m.group(1).lower()
        b64 = m.group(2).replace("\n", "")
        for key in ["org", "loo", "within", "delta", "meta"]:
            if key in tag:
                figs[key] = b64
                break
    return figs


def b64_to_img(b64: str):
    data = base64.b64decode(b64)
    return mpimg.imread(io.BytesIO(data), format="png")


def placeholder_img():
    arr = np.full((100, 200, 3), 0.92)
    return arr


# ── Slide chrome ──────────────────────────────────────────────────────────────
TOTAL_SLIDES = 10


def new_slide(title: str, n: int) -> plt.Figure:
    fig = plt.figure(figsize=(W, H), facecolor=WHITE)

    # Header bar
    ax_h = fig.add_axes([0, 0.875, 1, 0.125])
    ax_h.set_facecolor(NAVY)
    ax_h.set_xlim(0, 1); ax_h.set_ylim(0, 1); ax_h.axis("off")
    ax_h.text(0.025, 0.5, title, color=WHITE, fontsize=15, fontweight="bold", va="center")
    ax_h.text(0.975, 0.5, f"{n} / {TOTAL_SLIDES}", color=WHITE, fontsize=10,
              va="center", ha="right", alpha=0.65)

    # Orange accent stripe under header
    ax_s = fig.add_axes([0, 0.857, 1, 0.018])
    ax_s.set_facecolor(ORANGE); ax_s.axis("off")

    # Footer bar
    ax_f = fig.add_axes([0, 0, 1, 0.045])
    ax_f.set_facecolor(NAVY); ax_f.axis("off")
    ax_f.text(0.5, 0.5, "ProSurf  ·  Cross-Species OGT Analysis  ·  2026",
              color=WHITE, fontsize=8, va="center", ha="center", alpha=0.55)

    return fig


# ── Slide 1: Title ────────────────────────────────────────────────────────────
def slide_title(pdf: PdfPages):
    fig = plt.figure(figsize=(W, H), facecolor=NAVY)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_facecolor(NAVY); ax.axis("off")

    # Accent bar
    ax.axhline(0.44, xmin=0.07, xmax=0.93, color=ORANGE, linewidth=3)

    ax.text(0.5, 0.80, "Zwitterionic Surface Charge Mixing",
            color=WHITE, fontsize=26, fontweight="bold", ha="center", va="center",
            transform=ax.transAxes)
    ax.text(0.5, 0.68, "as a Predictor of Protein Thermal Stability",
            color=WHITE, fontsize=26, fontweight="bold", ha="center", va="center",
            transform=ax.transAxes)
    ax.text(0.5, 0.56, "A Cross-Species Analysis",
            color=ORANGE, fontsize=18, ha="center", va="center", transform=ax.transAxes)
    ax.text(0.5, 0.35, "13 organisms  ·  OGT 15–70°C  ·  34,231 AlphaFold2 protein structures",
            color="#aec6e8", fontsize=13, ha="center", va="center", transform=ax.transAxes)
    ax.text(0.5, 0.25, "Meltome Atlas (Jarzab et al., Nature Methods 2020)  |  ProSurf Z Score",
            color="#aec6e8", fontsize=11, ha="center", va="center", transform=ax.transAxes)
    ax.text(0.5, 0.14, "June 2026",
            color="#aec6e8", fontsize=11, ha="center", va="center", transform=ax.transAxes)

    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)


# ── Slide 2: Background ───────────────────────────────────────────────────────
def slide_background(pdf: PdfPages):
    fig = new_slide("Background & Motivation", 2)

    # ── Left: Z score schematic ──
    ax_l = fig.add_axes([0.03, 0.10, 0.44, 0.73])
    ax_l.set_facecolor(LGREY); ax_l.set_xlim(0, 10); ax_l.set_ylim(0, 10); ax_l.axis("off")

    ax_l.text(5, 9.4, "The ProSurf Zwitterionic Surface Score", fontsize=10.5,
              fontweight="bold", color=NAVY, ha="center")

    # Protein ellipse
    ell = Ellipse((5, 5.8), 7.2, 4.3, facecolor="#dce8f5", edgecolor=BLUE, linewidth=2)
    ax_l.add_patch(ell)
    angles = np.linspace(0, 2 * np.pi, 12, endpoint=False)
    for i, theta in enumerate(angles):
        x = 3.6 * np.cos(theta) + 5
        y = 2.15 * np.sin(theta) + 5.8
        ch, col = ("+", RED) if i % 2 == 0 else ("−", BLUE)
        ax_l.text(x, y, ch, color=col, fontsize=18, fontweight="bold", ha="center", va="center")

    ax_l.text(5, 5.8, "Z = B · D̂ · M", fontsize=15, fontweight="bold",
              color=DARK, ha="center", va="center")
    ax_l.text(5, 3.1, "Balance  ×  Density  ×  Mixing", fontsize=10, color="grey", ha="center")

    for row_y, lbl, desc in [
        (2.3, "B  Balance:", "min(n₊,n₋) / max(n₊,n₋)"),
        (1.55, "D̂  Density:", "charged residues / SASA  [normalised]"),
        (0.80, "M  Mixing:", "fraction of opposite-sign neighbours"),
    ]:
        ax_l.text(0.3, row_y, lbl, fontsize=9, fontweight="bold", color=DARK)
        ax_l.text(3.3, row_y, desc, fontsize=9, color="grey")

    # ── Right: Hypotheses ──
    ax_r = fig.add_axes([0.52, 0.10, 0.46, 0.73])
    ax_r.set_facecolor(WHITE); ax_r.set_xlim(0, 10); ax_r.set_ylim(0, 10); ax_r.axis("off")
    ax_r.text(5, 9.4, "Research Questions", fontsize=11, fontweight="bold",
              color=NAVY, ha="center")

    hyp_data = [
        (BLUE,   "H1", "Do organisms adapted to higher OGT have\nproteomes with higher mean Z?"),
        (GREEN,  "H2", "Within each organism, do thermally stable\nproteins have higher Z?"),
        (ORANGE, "H3", "Is any Z–Tm relationship independent\nof bulk charge composition?"),
    ]
    y0 = 8.2
    for col, tag, text in hyp_data:
        ax_r.add_patch(FancyBboxPatch((0.2, y0 - 0.45), 9.4, 1.3,
                                      boxstyle="round,pad=0.1",
                                      facecolor=col, alpha=0.10,
                                      edgecolor=col, linewidth=1.5))
        ax_r.text(0.7, y0 + 0.2, tag, fontsize=13, fontweight="bold", color=col, va="center")
        ax_r.text(2.2, y0 + 0.2, text, fontsize=9.5, color=DARK, va="center")
        y0 -= 2.1

    ax_r.text(5, 2.6, "Data sources", fontsize=10, fontweight="bold", color=NAVY, ha="center")
    for dy, line in enumerate([
        "Meltome Atlas — 13 organisms, per-protein Tm (TPP)",
        "AlphaFold2 v6 structures via EBI REST API",
        "34,231 proteins after quality filtering (≥150 residues)",
    ]):
        ax_r.text(5, 2.0 - dy * 0.6, line, fontsize=9, color="grey", ha="center")

    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)


# ── Slide 3: Dataset ──────────────────────────────────────────────────────────
def slide_dataset(pdf: PdfPages):
    fig = new_slide("Dataset: 13 Organisms · OGT 15–70°C · 34,231 Proteins", 3)
    ax = fig.add_axes([0.03, 0.07, 0.94, 0.77])
    ax.set_facecolor(WHITE); ax.set_xlim(0, 13); ax.set_ylim(0, 10); ax.axis("off")

    rows = [
        ("Oleispira antarctica",    15,  1261,  "Bacteria",  "#3498db"),
        ("Caenorhabditis elegans",  20,  3032,  "Eukaryota", "#27ae60"),
        ("Arabidopsis thaliana",    25,  2380,  "Eukaryota", "#27ae60"),
        ("Danio rerio",             28,  1786,  "Eukaryota", "#27ae60"),
        ("Drosophila melanogaster", 28,  1560,  "Eukaryota", "#27ae60"),
        ("Bacillus subtilis",       30,  1354,  "Bacteria",  "#3498db"),
        ("Saccharomyces cerevisiae",30,  1757,  "Eukaryota", "#27ae60"),
        ("Escherichia coli",        37,  1987,  "Bacteria",  "#3498db"),
        ("Homo sapiens",            37, 10695,  "Eukaryota", "#27ae60"),
        ("Mus musculus",            37,  6434,  "Eukaryota", "#27ae60"),
        ("Picrophilus torridus",    60,   767,  "Archaea",   "#e67e22"),
        ("Thermus thermophilus",    70,  1118,  "Bacteria",  "#3498db"),
    ]

    # Column headers
    for x, lbl in [(0.3, "Organism (italic)"), (6.5, "OGT (°C)"),
                    (7.9, "# Proteins"), (9.2, "Domain"), (10.6, "OGT range →")]:
        ax.text(x, 9.6, lbl, fontsize=9.5, fontweight="bold", color=NAVY)
    ax.axhline(9.2, xmin=0.02, xmax=0.98, color=NAVY, linewidth=0.8, alpha=0.4)

    max_bar = 2.5
    for i, (org, ogt, n_prot, domain, col) in enumerate(rows):
        y = 8.8 - i * 0.72
        ax.text(0.2, y, f"• {org}", fontsize=8.5, style="italic", color=DARK, va="center")
        ax.text(6.7, y, str(ogt), fontsize=9, color=col, fontweight="bold",
                ha="center", va="center")
        ax.text(8.1, y, f"{n_prot:,}", fontsize=9, color=DARK, ha="center", va="center")
        ax.text(9.4, y, domain, fontsize=8.5, color=col, va="center")
        bar_w = (ogt - 15) / (70 - 15) * max_bar
        ax.add_patch(plt.Rectangle((10.5, y - 0.25), bar_w, 0.5, color=col, alpha=0.75))

    ax.text(10.5, 9.6, "15°C", fontsize=7.5, color="grey")
    ax.text(13.0, 9.6, "70°C", fontsize=7.5, color="grey", ha="right")

    legend_handles = [
        mpatches.Patch(color="#3498db", label="Bacteria (4)"),
        mpatches.Patch(color="#27ae60", label="Eukaryota (7)"),
        mpatches.Patch(color="#e67e22", label="Archaea (1)"),
    ]
    ax.legend(handles=legend_handles, loc="lower right", fontsize=9, framealpha=0.85)

    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)


# ── Generic figure + bullet slide ─────────────────────────────────────────────
def slide_fig_bullets(pdf: PdfPages, n: int, title: str,
                      b64: str, bullets: list[str], key_col=GREEN):
    fig = new_slide(title, n)

    # Left: embedded analysis figure
    ax_img = fig.add_axes([0.01, 0.09, 0.54, 0.75])
    try:
        img = b64_to_img(b64) if b64 else placeholder_img()
    except Exception:
        img = placeholder_img()
    ax_img.imshow(img, aspect="auto")
    ax_img.axis("off")

    # Right: bullet points
    ax_r = fig.add_axes([0.57, 0.09, 0.42, 0.75])
    ax_r.set_facecolor(WHITE); ax_r.set_xlim(0, 10); ax_r.set_ylim(0, 10); ax_r.axis("off")

    y = 9.6
    for bullet in bullets:
        if bullet.startswith("##"):  # section heading
            ax_r.text(0, y, bullet[2:].strip(), fontsize=10.5, fontweight="bold",
                      color=NAVY, va="top")
            y -= 1.0
        elif bullet.startswith("!!"):  # key finding box
            txt = bullet[2:].strip()
            ax_r.add_patch(FancyBboxPatch((0, y - 0.7), 9.8, 0.85,
                                          boxstyle="round,pad=0.1",
                                          facecolor=key_col, alpha=0.15,
                                          edgecolor=key_col, linewidth=1.5))
            ax_r.text(0.3, y - 0.28, "▶  " + txt, fontsize=9.5,
                      color=key_col, fontweight="bold", va="center")
            y -= 1.1
        else:
            # Soft wrap at ~44 chars
            words = bullet.split()
            lines, cur = [], ""
            for w in words:
                if len(cur) + len(w) + 1 > 44:
                    lines.append(cur)
                    cur = w
                else:
                    cur = (cur + " " + w).strip()
            if cur:
                lines.append(cur)
            for j, line in enumerate(lines):
                pre = "•  " if j == 0 else "   "
                ax_r.text(0.3, y, pre + line, fontsize=9.5, color=DARK, va="top")
                y -= 0.56
            y -= 0.12

    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)


# ── Slide 9: Conclusions ──────────────────────────────────────────────────────
def slide_conclusions(pdf: PdfPages):
    fig = new_slide("Conclusions", 9)
    ax = fig.add_axes([0.03, 0.07, 0.94, 0.77])
    ax.set_facecolor(WHITE); ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

    items = [
        (BLUE,   "H1  Cross-species OGT trend", "NOT SUPPORTED",
         "Organism-level ρ(Z, OGT) = +0.21, permutation p = 0.50 (N = 13 organisms).\n"
         "Signal collapses to ρ = 0 when either thermophile is removed (LOO test).\n"
         "After composition control: partial ρ = +0.006, p = 0.30."),
        (GREEN,  "H2  Within-organism Tm association", "SUPPORTED",
         "Pooled ρ(Z, Tm) = +0.079, p = 1.26×10⁻⁴⁸  (N = 34,231 proteins).\n"
         "Strongest in bacteria: B. subtilis +0.149, E. coli +0.109.\n"
         "Absent in Picrophilus (+0.031, n.s.) and Drosophila (−0.026, n.s.)."),
        (ORANGE, "H3  Composition independence", "SUPPORTED",
         "Z adds ΔR² = +0.0038 (p = 3.67×10⁻³⁰) beyond bulk charge features.\n"
         "Coefficient: +0.688 °C per σ_Z (semi-standardised).\n"
         "Z–Tm coupling strength does not scale with OGT (meta-regression p = 0.71)."),
    ]

    y = 9.1
    for col, heading, verdict, detail in items:
        v_col = RED if "NOT" in verdict else GREEN
        ax.add_patch(FancyBboxPatch((0.1, y - 2.0), 9.8, 2.3,
                                    boxstyle="round,pad=0.1",
                                    facecolor=col, alpha=0.07,
                                    edgecolor=col, linewidth=1.2))
        ax.text(0.4, y - 0.3, heading, fontsize=11, fontweight="bold", color=col, va="center")
        ax.text(9.7, y - 0.3, verdict, fontsize=10, fontweight="bold",
                color=v_col, va="center", ha="right")
        for k, line in enumerate(detail.split("\n")):
            ax.text(0.6, y - 0.85 - k * 0.55, "•  " + line.strip(),
                    fontsize=9, color=DARK, va="top")
        y -= 3.05

    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)


# ── Slide 10: Take-home ───────────────────────────────────────────────────────
def slide_takehome(pdf: PdfPages):
    fig = new_slide("Take-Home Message", 10)
    ax = fig.add_axes([0.03, 0.07, 0.94, 0.77])
    ax.set_facecolor(WHITE); ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

    panels = [
        (GREEN, "#e8f5ee",
         "Within proteomes  ✓",
         "Proteins with better-mixed surface charges\n"
         "(checkerboard-like Z) are measurably more\n"
         "thermally stable.\n\n"
         "~0.7 °C per σ_Z,  independent of total\n"
         "charge count.  (p = 3.67×10⁻³⁰)"),
        (RED, "#fde8e8",
         "Across organisms  ✗",
         "Thermophilic proteomes do NOT show\n"
         "higher Z than mesophiles once bulk\n"
         "charge density is accounted for.\n\n"
         "OGT adaptation uses more charges,\n"
         "not better charge mixing."),
    ]

    for i, (col, bg, ttl, body) in enumerate(panels):
        x0 = 0.3 + i * 5.0
        ax.add_patch(FancyBboxPatch((x0, 2.5), 4.5, 6.5,
                                    boxstyle="round,pad=0.15",
                                    facecolor=bg, edgecolor=col, linewidth=2))
        ax.text(x0 + 2.25, 8.6, ttl, fontsize=12.5, fontweight="bold",
                color=col, ha="center")
        ax.axhline(8.1, xmin=(x0 + 0.1) / 10, xmax=(x0 + 4.4) / 10,
                   color=col, linewidth=1, alpha=0.4)
        for j, line in enumerate(body.split("\n")):
            ax.text(x0 + 0.25, 7.55 - j * 0.63, line, fontsize=9.5, color=DARK, va="top")

    ax.text(5, 2.0,
            "Engineering implication: improving surface charge mixing (Z) is a real,"
            " composition-independent lever for\n"
            "single-protein stabilisation — particularly for bacterial enzyme targets."
            "  Combine with core-packing and disulfide strategies.",
            fontsize=9.5, color=NAVY, ha="center", va="top",
            style="italic", linespacing=1.5)

    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    print(f"Extracting figures from {SOURCE_REPORT} ...")
    figs = extract_figures(SOURCE_REPORT)
    found = list(figs.keys())
    print(f"  Found: {found}")
    missing = [k for k in ["org", "loo", "within", "delta", "meta"] if k not in figs]
    if missing:
        print(f"  Warning: missing figure keys: {missing} — placeholders will be used")

    print(f"Generating {TOTAL_SLIDES}-slide deck ...")

    with PdfPages(OUT_PDF) as pdf:
        d = pdf.infodict()
        d["Title"]    = "Zwitterionic Surface Score and Protein Thermal Stability"
        d["Subject"]  = "Cross-Species OGT Analysis"
        d["Keywords"] = "ProSurf, Z score, thermal stability, OGT, Meltome Atlas"

        slide_title(pdf)                                                        # 1
        slide_background(pdf)                                                   # 2
        slide_dataset(pdf)                                                      # 3

        slide_fig_bullets(pdf, 4,                                               # 4
            "Result 1 — Organism-level Z vs OGT",
            figs.get("org", ""),
            [
                "## Spearman ρ(mean Z, OGT) = +0.21,  perm-p = 0.50",
                "N = 13 organisms, OGT 15–70°C",
                "Permutation test with 10,000 shuffles",
                "Consistent across z_mean / z_max / z_frac (ρ = +0.15–0.21)",
                "Stable across MIN_PER_ORG cutoffs 30–200",
                "10 effective phylogenetic groups — inherently low power",
                "!!H1 NOT supported: no significant cross-species OGT trend",
            ],
            key_col=RED,
        )

        slide_fig_bullets(pdf, 5,                                               # 5
            "Result 2 — Leave-One-Out Robustness",
            figs.get("loo", ""),
            [
                "## LOO ρ collapses to zero without thermophiles",
                "Drop Picrophilus (60°C)  →  ρ = 0.000",
                "Drop Thermus (70°C)       →  ρ = 0.000",
                "All 11 mesophile-only subsets: ρ = +0.175 to +0.329",
                "Trend is a two-point extrapolation, not a broad signal",
                "Drop E. coli alone reduces ρ by 62% (anomalously low Z)",
                "!!Cross-species trend depends entirely on 2 thermophiles",
            ],
            key_col=RED,
        )

        slide_fig_bullets(pdf, 6,                                               # 6
            "Result 3 — Within-Organism Z–Tm Association",
            figs.get("within", ""),
            [
                "## Pooled ρ = +0.079,  p = 1.26×10⁻⁴⁸  (N = 34,231)",
                "11 of 13 organisms show positive trend",
                "Bacteria strongest: B. subtilis +0.149, E. coli +0.109",
                "Vertebrates moderate: H. sapiens +0.086, M. musculus +0.081",
                "Drosophila: −0.026 (n.s.) — only negative outlier",
                "Psychrophile Oleispira & acidophile Picrophilus: n.s.",
                "!!H2 SUPPORTED: within-proteome Z predicts Tm",
            ],
            key_col=GREEN,
        )

        slide_fig_bullets(pdf, 7,                                               # 7
            "Result 4 — Composition Control",
            figs.get("delta", ""),
            [
                "## Z adds beyond bulk charge features",
                "Partial ρ(Z, OGT | composition) = +0.006,  p = 0.30",
                "  → Z carries no independent OGT signal",
                "Nested Tm model:  ΔR² = +0.0038,  p = 3.67×10⁻³⁰",
                "Semi-standardised coef: +0.688 °C / σ_Z",
                "Effect size small (0.38% variance) — robust at N = 34k",
                "!!H3 SUPPORTED: Z–Tm is composition-independent",
            ],
            key_col=GREEN,
        )

        slide_fig_bullets(pdf, 8,                                               # 8
            "Result 5 — Meta-regression: Domain-of-Life Specificity",
            figs.get("meta", ""),
            [
                "## Z–Tm coupling does not scale with OGT",
                "Meta-regression ρ(ρ_within, OGT) = +0.117,  perm-p = 0.71",
                "Domain-level gradient (Z as stability predictor):",
                "  Bacteria  :  mean ρ = +0.096  (n = 4 organisms)",
                "  Eukaryota :  mean ρ = +0.044  (n = 8 organisms)",
                "  Archaea   :  mean ρ = +0.031  (n = 1 organism)",
                "!!Z is most predictive for bacterial protein targets",
            ],
            key_col=BLUE,
        )

        slide_conclusions(pdf)                                                   # 9
        slide_takehome(pdf)                                                      # 10

    size_kb = OUT_PDF.stat().st_size // 1024
    print(f"Wrote {OUT_PDF}  ({size_kb} KB,  {TOTAL_SLIDES} slides)")


if __name__ == "__main__":
    main()
