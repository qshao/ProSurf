#!/usr/bin/env python3
"""
Generate a professional summary slide deck PDF for the ProSurf cross-species OGT analysis.

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

# Slide canvas: 16:9
W, H = 13.33, 7.5

# ── Palette ───────────────────────────────────────────────────────────────────
NAVY    = "#15304F"
STEEL   = "#2472A4"
TEAL    = "#0E8C8C"
PINE    = "#18795A"
CRIMSON = "#9B2335"
AMBER   = "#C96B10"
GREY    = "#5A6A7A"
LGREY   = "#F5F7FA"
MID     = "#DDE3EA"
WHITE   = "#FFFFFF"
INK     = "#1A2030"

TOTAL = 10


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
    return mpimg.imread(io.BytesIO(base64.b64decode(b64)), format="png")


def fit_image(fig_w, fig_h, rl, rb, rw, rh, img_ar, pad=0.01):
    """Return axes [l, b, w, h] in figure fractions that fit img_ar inside region without distortion."""
    aw = (rw - 2 * pad) * fig_w
    ah = (rh - 2 * pad) * fig_h
    if aw / ah > img_ar:
        ih = ah; iw = ih * img_ar
    else:
        iw = aw; ih = iw / img_ar
    l = rl + pad + (aw - iw) / (2 * fig_w)
    b = rb + pad + (ah - ih) / (2 * fig_h)
    return [l, b, iw / fig_w, ih / fig_h]


# ── Chrome helpers ────────────────────────────────────────────────────────────
def chrome(fig: plt.Figure, title: str, n: int):
    """Draw header bar, accent stripe and footer on fig."""
    # Header
    ax = fig.add_axes([0, 0.905, 1, 0.095])
    ax.set_facecolor(NAVY); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.022, 0.48, title, color=WHITE, fontsize=14, fontweight="bold", va="center")
    ax.text(0.978, 0.48, f"{n} / {TOTAL}", color=WHITE, fontsize=9,
            va="center", ha="right", alpha=0.6)
    # Teal accent stripe
    ax2 = fig.add_axes([0, 0.893, 1, 0.012])
    ax2.set_facecolor(TEAL); ax2.axis("off")
    # Footer
    ax3 = fig.add_axes([0, 0, 1, 0.038])
    ax3.set_facecolor(NAVY); ax3.axis("off")
    ax3.text(0.5, 0.5,
             "ProSurf  ·  Cross-Species Zwitterionic Surface Analysis  ·  2026",
             color=WHITE, fontsize=7.5, va="center", ha="center", alpha=0.5)


def new_slide(title: str, n: int) -> plt.Figure:
    fig = plt.figure(figsize=(W, H), facecolor=WHITE)
    chrome(fig, title, n)
    return fig


def content_region():
    """Standard content area: [left, bottom, width, height] in figure fractions."""
    return 0.015, 0.055, 0.970, 0.830


# ── Metric card helper ────────────────────────────────────────────────────────
def metric_card(ax, x, y, w, h, label, value, col):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05",
                                 facecolor=col, alpha=0.10, edgecolor=col, linewidth=1.2))
    ax.text(x + w / 2, y + h * 0.68, value, fontsize=11.5, fontweight="bold",
            color=col, ha="center", va="center")
    ax.text(x + w / 2, y + h * 0.22, label, fontsize=7.5,
            color=GREY, ha="center", va="center")


def verdict_box(ax, x, y, w, text, col):
    h = 0.9
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08",
                                 facecolor=col, alpha=0.13, edgecolor=col, linewidth=1.8))
    ax.text(x + w / 2, y + h / 2, text, fontsize=9.5, fontweight="bold",
            color=col, ha="center", va="center", linespacing=1.4)


def bullets(ax, x, y0, items, fontsize=9.5, dy=0.65):
    y = y0
    for item in items:
        ax.text(x, y, "▸  " + item, fontsize=fontsize, color=INK, va="top", linespacing=1.3)
        y -= dy + item.count("\n") * 0.45
    return y


# ── Slide 1: Title ────────────────────────────────────────────────────────────
def slide_title(pdf: PdfPages):
    fig = plt.figure(figsize=(W, H), facecolor=NAVY)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_facecolor(NAVY); ax.axis("off")

    # Subtle geometric accent — diagonal stripe band
    for xi in np.linspace(-1, 2, 30):
        ax.plot([xi, xi + 0.6], [0, 1], color=TEAL, alpha=0.04,
                linewidth=8, transform=ax.transAxes)

    # Teal horizontal rule
    ax.plot([0.07, 0.93], [0.40, 0.40], color=TEAL, linewidth=2.5,
            transform=ax.transAxes)

    ax.text(0.5, 0.82, "Zwitterionic Surface Charge Mixing",
            color=WHITE, fontsize=27, fontweight="bold",
            ha="center", va="center", transform=ax.transAxes)
    ax.text(0.5, 0.70, "as a Predictor of Protein Thermal Stability",
            color=WHITE, fontsize=27, fontweight="bold",
            ha="center", va="center", transform=ax.transAxes)
    ax.text(0.5, 0.58, "A Cross-Species Analysis",
            color=TEAL, fontsize=16, ha="center", va="center", transform=ax.transAxes)

    ax.text(0.5, 0.31,
            "13 organisms  ·  OGT 15–70°C  ·  34,231 AlphaFold2 protein structures",
            color="#8BAEC8", fontsize=12.5, ha="center", va="center",
            transform=ax.transAxes)
    ax.text(0.5, 0.21,
            "Meltome Atlas (Jarzab et al., Nature Methods 2020)",
            color="#8BAEC8", fontsize=11, ha="center", va="center",
            transform=ax.transAxes)
    ax.text(0.5, 0.10, "June 2026",
            color="#8BAEC8", fontsize=11, ha="center", va="center",
            transform=ax.transAxes)

    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)


# ── Slide 2: Background ───────────────────────────────────────────────────────
def slide_background(pdf: PdfPages):
    fig = new_slide("Background & Motivation", 2)

    # Left panel – Z schematic
    ax_l = fig.add_axes([0.02, 0.06, 0.43, 0.82])
    ax_l.set_facecolor(LGREY); ax_l.set_xlim(0, 10); ax_l.set_ylim(0, 11)
    ax_l.axis("off")

    ax_l.text(5, 10.4, "The ProSurf Zwitterionic Surface Score",
              fontsize=10, fontweight="bold", color=NAVY, ha="center")

    # Protein ellipse
    ell = Ellipse((5, 6.2), 7.0, 4.0,
                  facecolor="#D6E8F5", edgecolor=STEEL, linewidth=1.8)
    ax_l.add_patch(ell)
    for i, theta in enumerate(np.linspace(0, 2 * np.pi, 12, endpoint=False)):
        x = 3.5 * np.cos(theta) + 5
        y = 2.0 * np.sin(theta) + 6.2
        ch, col = ("+", "#C0392B") if i % 2 == 0 else ("−", STEEL)
        ax_l.text(x, y, ch, color=col, fontsize=17, fontweight="bold",
                  ha="center", va="center")

    ax_l.text(5, 6.2, "Z = B · D̂ · M", fontsize=14, fontweight="bold",
              color=INK, ha="center", va="center")
    ax_l.text(5, 3.55, "Balance × Density × Mixing", fontsize=9.5,
              color=GREY, ha="center")

    sep_y = 3.0
    ax_l.axhline(sep_y, xmin=0.05, xmax=0.95, color=MID, linewidth=1)

    for row_y, tag, lbl in [
        (2.4, "B", "Balance   min(n₊,n₋) / max(n₊,n₋)"),
        (1.7, "D̂", "Density   charged residues / SASA  [norm.]"),
        (1.0, "M", "Mixing    fraction of opposite-sign neighbours"),
    ]:
        ax_l.text(0.5, row_y, tag, fontsize=10, fontweight="bold",
                  color=TEAL, va="center")
        ax_l.text(1.5, row_y, lbl, fontsize=8.5, color=GREY, va="center")

    # Right panel – Hypotheses
    ax_r = fig.add_axes([0.48, 0.06, 0.50, 0.82])
    ax_r.set_facecolor(WHITE); ax_r.set_xlim(0, 10); ax_r.set_ylim(0, 11)
    ax_r.axis("off")

    ax_r.text(5, 10.4, "Research Questions", fontsize=11,
              fontweight="bold", color=NAVY, ha="center")

    for y0, col, tag, q in [
        (8.5, STEEL,   "H1", "Do organisms adapted to higher OGT have\nproteomes with higher mean Z?"),
        (5.9, PINE,    "H2", "Within each organism, do thermally stable\nproteins carry higher Z scores?"),
        (3.3, AMBER,   "H3", "Is any Z–Tm relationship independent\nof bulk charge composition?"),
    ]:
        ax_r.add_patch(FancyBboxPatch((0.1, y0 - 0.5), 9.6, 2.0,
                                       boxstyle="round,pad=0.1",
                                       facecolor=col, alpha=0.08,
                                       edgecolor=col, linewidth=1.3))
        ax_r.text(0.65, y0 + 0.55, tag, fontsize=13, fontweight="bold",
                  color=col, va="center")
        ax_r.text(2.2, y0 + 0.55, q, fontsize=9.5, color=INK, va="center",
                  linespacing=1.5)

    ax_r.text(5, 2.4, "Data", fontsize=9, fontweight="bold", color=NAVY, ha="center")
    for i, line in enumerate([
        "Meltome Atlas — 13 organisms, per-protein Tm (thermal proteome profiling)",
        "AlphaFold2 v6 structures via EBI REST API",
        "34,231 proteins after quality filtering  (≥ 150 residues, ≥ 30 per organism)",
    ]):
        ax_r.text(5, 1.85 - i * 0.55, line, fontsize=8.5,
                  color=GREY, ha="center")

    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)


# ── Slide 3: Dataset ──────────────────────────────────────────────────────────
def slide_dataset(pdf: PdfPages):
    fig = new_slide("Dataset Overview", 3)

    ax = fig.add_axes([0.02, 0.06, 0.96, 0.82])
    ax.set_facecolor(WHITE); ax.set_xlim(0, 14); ax.set_ylim(0, 11); ax.axis("off")

    rows = [
        ("Oleispira antarctica",     15,  1261, "Bacteria",  STEEL),
        ("Caenorhabditis elegans",   20,  3032, "Eukaryota", PINE),
        ("Arabidopsis thaliana",     25,  2380, "Eukaryota", PINE),
        ("Danio rerio",              28,  1786, "Eukaryota", PINE),
        ("Drosophila melanogaster",  28,  1560, "Eukaryota", PINE),
        ("Bacillus subtilis",        30,  1354, "Bacteria",  STEEL),
        ("Saccharomyces cerevisiae", 30,  1757, "Eukaryota", PINE),
        ("Escherichia coli",         37,  1987, "Bacteria",  STEEL),
        ("Homo sapiens",             37, 10695, "Eukaryota", PINE),
        ("Mus musculus",             37,  6434, "Eukaryota", PINE),
        ("Picrophilus torridus",     60,   767, "Archaea",   AMBER),
        ("Thermus thermophilus",     70,  1118, "Bacteria",  STEEL),
    ]

    # Header row
    for x, lbl, ha in [
        (0.2, "Organism", "left"),
        (6.6, "OGT (°C)", "center"),
        (7.9, "Proteins", "center"),
        (9.2, "Domain", "left"),
        (10.8, "← 15°C", "left"),
        (13.7, "70°C →", "right"),
    ]:
        ax.text(x, 10.5, lbl, fontsize=9, fontweight="bold", color=NAVY, ha=ha)
    ax.axhline(10.1, xmin=0.01, xmax=0.99, color=MID, linewidth=1.2)

    max_bar = 2.8
    for i, (org, ogt, n_prot, domain, col) in enumerate(rows):
        y = 9.55 - i * 0.76
        bg = LGREY if i % 2 == 0 else WHITE
        ax.add_patch(plt.Rectangle((0.05, y - 0.34), 13.9, 0.70,
                                    facecolor=bg, zorder=0))
        ax.text(0.25, y, f"• {org}", fontsize=8.5, style="italic",
                color=INK, va="center")
        ax.text(6.6, y, str(ogt), fontsize=9, color=col,
                fontweight="bold", ha="center", va="center")
        ax.text(7.9, y, f"{n_prot:,}", fontsize=9, color=INK,
                ha="center", va="center")
        ax.text(9.2, y, domain, fontsize=8.5, color=col, va="center")
        bar_w = (ogt - 15) / (70 - 15) * max_bar
        ax.add_patch(plt.Rectangle((10.8, y - 0.22), bar_w, 0.44,
                                    facecolor=col, alpha=0.70, zorder=1))

    legend_handles = [
        mpatches.Patch(color=STEEL, label="Bacteria (4)"),
        mpatches.Patch(color=PINE,  label="Eukaryota (7)"),
        mpatches.Patch(color=AMBER, label="Archaea (1)"),
    ]
    ax.legend(handles=legend_handles, loc="lower right",
              fontsize=8.5, framealpha=0.9, edgecolor=MID)

    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)


# ── Generic figure + analysis slide ──────────────────────────────────────────
def slide_result(pdf: PdfPages, n: int, title: str, img_ar: float, b64: str,
                 metrics: list[tuple],   # [(label, value, col), ...]
                 bullet_items: list[str],
                 verdict_text: str, verdict_col: str):
    fig = new_slide(title, n)

    # ── Left: analysis figure ────────────────────────────────────────────────
    # Light card background
    ax_card = fig.add_axes([0.015, 0.055, 0.550, 0.840])
    ax_card.set_facecolor(LGREY); ax_card.axis("off")

    rect = fit_image(W, H, 0.015, 0.055, 0.550, 0.840, img_ar, pad=0.012)
    ax_img = fig.add_axes(rect)
    try:
        img = b64_to_img(b64)
        ax_img.imshow(img)
    except Exception:
        ax_img.set_facecolor(MID)
    ax_img.axis("off")

    # ── Right: metrics + bullets ─────────────────────────────────────────────
    ax_r = fig.add_axes([0.585, 0.055, 0.400, 0.840])
    ax_r.set_facecolor(WHITE); ax_r.set_xlim(0, 10); ax_r.set_ylim(0, 11)
    ax_r.axis("off")

    # Metric cards (2–3 small KPI boxes at top)
    card_w = 9.6 / len(metrics)
    for j, (lbl, val, col) in enumerate(metrics):
        metric_card(ax_r, 0.2 + j * card_w, 8.8, card_w - 0.3, 1.5,
                    lbl, val, col)

    # Separator
    ax_r.axhline(8.4, xmin=0, xmax=1, color=MID, linewidth=0.8)

    # Bullet points
    y = 7.9
    for item in bullet_items:
        # word-wrap at ~42 chars
        words = item.split()
        lines, cur = [], ""
        for w in words:
            if len(cur) + len(w) + 1 > 42:
                lines.append(cur); cur = w
            else:
                cur = (cur + " " + w).strip()
        if cur:
            lines.append(cur)
        for k, line in enumerate(lines):
            pre = "▸  " if k == 0 else "    "
            ax_r.text(0.2, y, pre + line, fontsize=9.5, color=INK, va="top")
            y -= 0.56
        y -= 0.10

    # Verdict box at bottom
    verdict_box(ax_r, 0.1, 0.4, 9.6, verdict_text, verdict_col)

    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)


# ── Slide 2 (Background) already defined above ────────────────────────────────

# ── Slide 9: Conclusions ──────────────────────────────────────────────────────
def slide_conclusions(pdf: PdfPages):
    fig = new_slide("Conclusions", 9)
    ax = fig.add_axes([0.015, 0.055, 0.970, 0.840])
    ax.set_facecolor(WHITE); ax.set_xlim(0, 10); ax.set_ylim(0, 11); ax.axis("off")

    verdicts = [
        (STEEL,   "H1  Cross-species OGT trend",     "✗  NOT SUPPORTED",  CRIMSON,
         ["Organism-level ρ(Z, OGT) = +0.21,  permutation p = 0.50",
          "Signal collapses to ρ = 0 when either thermophile is removed (LOO)",
          "After composition control: partial ρ = +0.006, p = 0.30"]),
        (PINE,    "H2  Within-organism Tm association", "✓  SUPPORTED",    PINE,
         ["Pooled ρ(Z, Tm) = +0.079, p = 1.26×10⁻⁴⁸  (N = 34,231 proteins)",
          "11 / 13 organisms show positive trend",
          "Strongest in bacteria: B. subtilis +0.149, E. coli +0.109"]),
        (AMBER,   "H3  Composition independence",     "✓  SUPPORTED",    PINE,
         ["Z adds ΔR² = +0.0038 (F-test p = 3.67×10⁻³⁰) beyond bulk charge",
          "Semi-standardised coefficient: +0.688 °C per σ_Z",
          "Z–Tm coupling does not scale with OGT (meta-regression p = 0.71)"]),
    ]

    y0 = 10.4
    for col, heading, verdict, v_col, detail in verdicts:
        ax.add_patch(FancyBboxPatch((0.1, y0 - 2.6), 9.7, 2.8,
                                    boxstyle="round,pad=0.12",
                                    facecolor=col, alpha=0.06,
                                    edgecolor=col, linewidth=1.4))
        ax.text(0.45, y0 - 0.5, heading, fontsize=11, fontweight="bold",
                color=col, va="center")
        ax.text(9.65, y0 - 0.5, verdict, fontsize=10, fontweight="bold",
                color=v_col, va="center", ha="right")
        ax.axhline(y0 - 0.95, xmin=0.03, xmax=0.97,
                   color=col, linewidth=0.7, alpha=0.4)
        for k, line in enumerate(detail):
            ax.text(0.6, y0 - 1.45 - k * 0.55, "•  " + line,
                    fontsize=9, color=INK, va="top")
        y0 -= 3.42

    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)


# ── Slide 10: Take-home ───────────────────────────────────────────────────────
def slide_takehome(pdf: PdfPages):
    fig = new_slide("Take-Home Message", 10)
    ax = fig.add_axes([0.015, 0.055, 0.970, 0.840])
    ax.set_facecolor(WHITE); ax.set_xlim(0, 10); ax.set_ylim(0, 11); ax.axis("off")

    panels = [
        (PINE,    "#E6F4ED", "Within individual proteomes  ✓",
         ["Proteins with more evenly mixed surface charges",
          "(higher Z) are measurably more thermally stable.",
          "",
          "Effect:  +0.7 °C per σ_Z",
          "p = 3.67×10⁻³⁰  (N = 34,231 proteins)",
          "Independent of bulk charge count",
          "",
          "Strongest in bacteria  (B. subtilis, E. coli)"]),
        (CRIMSON, "#F9EAEA", "Across organisms  ✗",
         ["Thermophilic proteomes do NOT show higher Z",
          "than mesophiles once bulk charge density",
          "is accounted for.",
          "",
          "Adaptation to high OGT → more charges,",
          "not better spatial charge mixing.",
          "",
          "Signal collapses without 2 thermophiles"]),
    ]

    for i, (col, bg, ttl, body) in enumerate(panels):
        x0 = 0.25 + i * 4.9
        ax.add_patch(FancyBboxPatch((x0, 1.6), 4.5, 8.3,
                                    boxstyle="round,pad=0.15",
                                    facecolor=bg, edgecolor=col, linewidth=2))
        ax.text(x0 + 2.25, 9.55, ttl, fontsize=11.5, fontweight="bold",
                color=col, ha="center")
        ax.axhline(9.1, xmin=(x0 + 0.1) / 10, xmax=(x0 + 4.4) / 10,
                   color=col, linewidth=1.2, alpha=0.5)
        y = 8.6
        for line in body:
            if line == "":
                y -= 0.3
            else:
                ax.text(x0 + 0.35, y, line, fontsize=9.5,
                        color=INK, va="top")
                y -= 0.59

    # Bottom engineering implication
    ax.add_patch(FancyBboxPatch((0.1, 0.1), 9.7, 1.2,
                                boxstyle="round,pad=0.1",
                                facecolor=TEAL, alpha=0.10,
                                edgecolor=TEAL, linewidth=1.5))
    ax.text(5.0, 0.72,
            "Engineering implication: improving surface charge mixing (Z) is a real,"
            " composition-independent lever for single-protein\n"
            "stabilisation — best validated for bacterial enzymes. Use alongside core-packing"
            " and disulfide strategies.",
            fontsize=9, color=NAVY, ha="center", va="center",
            linespacing=1.5, style="italic")

    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    print(f"Extracting figures from {SOURCE_REPORT.name} ...")
    figs = extract_figures(SOURCE_REPORT)
    print(f"  Found: {list(figs.keys())}")
    missing = [k for k in ["org", "loo", "within", "delta", "meta"] if k not in figs]
    if missing:
        print(f"  Warning — missing: {missing}")

    # Precompute aspect ratios (W/H) for each analysis figure
    ar = {}
    for k, b64 in figs.items():
        img = b64_to_img(b64)
        ar[k] = img.shape[1] / img.shape[0]
    print(f"  Aspect ratios: {', '.join(f'{k}={v:.2f}' for k,v in ar.items())}")

    print(f"Generating {TOTAL}-slide deck ...")

    with PdfPages(OUT_PDF) as pdf:
        d = pdf.infodict()
        d["Title"]    = "Zwitterionic Surface Score and Protein Thermal Stability"
        d["Subject"]  = "Cross-Species OGT Analysis · ProSurf"
        d["Keywords"] = "ProSurf, Z score, thermal stability, OGT, Meltome Atlas, AlphaFold2"

        # 1 — Title
        slide_title(pdf)

        # 2 — Background
        slide_background(pdf)

        # 3 — Dataset
        slide_dataset(pdf)

        # 4 — H1: organism-level
        slide_result(pdf, 4,
            "Result 1 — Organism-level Z vs Optimal Growth Temperature",
            ar.get("org", 1.32), figs.get("org", ""),
            metrics=[
                ("Spearman ρ",   "+0.21",  STEEL),
                ("Permutation p", "0.50",  CRIMSON),
                ("N organisms",  "13",     GREY),
            ],
            bullet_items=[
                "Positive but non-significant trend across OGT 15–70°C",
                "Consistent across all 3 Z metrics (ρ = +0.15 to +0.21)",
                "Stable across MIN_PER_ORG cutoffs 30–200",
                "~10 effective phylogenetic groups — low statistical power",
            ],
            verdict_text="H1 NOT SUPPORTED — no significant cross-species OGT trend",
            verdict_col=CRIMSON,
        )

        # 5 — LOO
        slide_result(pdf, 5,
            "Result 2 — Leave-One-Out Robustness of the OGT Trend",
            ar.get("loo", 1.28), figs.get("loo", ""),
            metrics=[
                ("ρ without Picrophilus", "0.000", CRIMSON),
                ("ρ without Thermus",     "0.000", CRIMSON),
                ("Mesophile-only range",  "+0.18–+0.33", GREY),
            ],
            bullet_items=[
                "Removing either thermophile collapses ρ exactly to zero",
                "Trend is a two-point extrapolation, not a broad signal",
                "Dropping E. coli alone reduces ρ by 62% (anomalously low Z)",
                "Dropping Homo sapiens gives highest LOO ρ = +0.329",
            ],
            verdict_text="Cross-species trend depends entirely on 2 thermophilic organisms",
            verdict_col=CRIMSON,
        )

        # 6 — Within-organism
        slide_result(pdf, 6,
            "Result 3 — Within-Organism Z–Tm Association",
            ar.get("within", 1.23), figs.get("within", ""),
            metrics=[
                ("Pooled ρ",    "+0.079", PINE),
                ("p-value",     "1.3×10⁻⁴⁸", PINE),
                ("N proteins",  "34,231", GREY),
            ],
            bullet_items=[
                "11 of 13 organisms show a positive Z–Tm relationship",
                "Bacteria strongest: B. subtilis +0.149, E. coli +0.109",
                "Vertebrates moderate: H. sapiens +0.086, M. musculus +0.081",
                "Psychrophile Oleispira (+0.038, n.s.) and acidophile Picrophilus (+0.031, n.s.) absent",
                "Drosophila (−0.026, n.s.) — sole negative outlier",
            ],
            verdict_text="H2 SUPPORTED — within-proteome Z robustly predicts Tm",
            verdict_col=PINE,
        )

        # 7 — Composition control
        slide_result(pdf, 7,
            "Result 4 — Composition Control (Partial Correlation & Nested OLS)",
            ar.get("delta", 1.05), figs.get("delta", ""),
            metrics=[
                ("ΔR²", "+0.0038", PINE),
                ("F-test p", "3.7×10⁻³⁰", PINE),
                ("Coef (+0.688 °C / σ_Z)", "composition-independent", TEAL),
            ],
            bullet_items=[
                "Partial ρ(Z, OGT | composition) = +0.006, p = 0.30 — OGT signal vanishes",
                "Z adds significant Tm prediction beyond charged_frac, net_charge, density",
                "Effect size small (0.38% variance) but robust at N = 34k",
                "+0.688 °C per σ_Z: modest but real, orthogonal to charge count",
            ],
            verdict_text="H3 SUPPORTED — Z–Tm association is composition-independent",
            verdict_col=PINE,
        )

        # 8 — Meta-regression
        slide_result(pdf, 8,
            "Result 5 — Meta-regression: Does Z–Tm Coupling Scale with OGT?",
            ar.get("meta", 1.34), figs.get("meta", ""),
            metrics=[
                ("ρ_meta",       "+0.117", GREY),
                ("Perm-p",       "0.71",   GREY),
                ("Domain gradient", "Bacteria > Eukaryota > Archaea", TEAL),
            ],
            bullet_items=[
                "Z–Tm coupling strength does NOT increase with OGT (p = 0.71)",
                "Bacteria: mean within-org ρ = +0.096  (n = 4 organisms)",
                "Eukaryota: mean within-org ρ = +0.044  (n = 8 organisms)",
                "Archaea: mean within-org ρ = +0.031  (n = 1 organism)",
                "Z is most predictive for bacterial enzyme targets",
            ],
            verdict_text="Z–Tm coupling is domain-specific, not temperature-tuned",
            verdict_col=TEAL,
        )

        # 9 — Conclusions
        slide_conclusions(pdf)

        # 10 — Take-home
        slide_takehome(pdf)

    size_kb = OUT_PDF.stat().st_size // 1024
    print(f"Wrote {OUT_PDF}  ({size_kb} KB,  {TOTAL} slides)")


if __name__ == "__main__":
    main()
