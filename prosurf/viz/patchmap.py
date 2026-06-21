from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def render_patch_map(arr, scores, out_png):
    out_png = Path(out_png)
    pts = np.array([s.xyz for s in scores]) if scores else np.zeros((1, 3))
    vals = np.array([s.z for s in scores]) if scores else np.array([0.0])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    p = ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], c=vals, cmap="viridis")
    fig.colorbar(p, label="Z (zwitterionic)")
    fig.savefig(out_png, dpi=120)
    plt.close(fig)
    return out_png
