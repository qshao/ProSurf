# scripts/run_validation.py
import yaml
from pathlib import Path
from prosurf.config import MetricConfig, PathsConfig
from prosurf.io.fetch import fetch_af2
from prosurf.pipeline import analyze_structure
from prosurf.validate.controls import auroc
from prosurf.validate.robustness import sweep_parameter, ranking_stability, engine_agreement

def main():
    cfg = MetricConfig(); paths = PathsConfig()
    cs = yaml.safe_load(Path("configs/control_set.yaml").read_text())
    def scores(accs):
        out = []
        for a in accs:
            p = fetch_af2(a, paths.data_dir)
            _, ps = analyze_structure(p, a, cfg)
            out.append(ps.z_max)
        return out
    pos, neg = scores(cs["positives"]), scores(cs["negatives"])
    au = auroc(pos, neg)
    structs = [(fetch_af2(a, paths.data_dir), a) for a in cs["positives"] + cs["negatives"]]
    lines = [f"# Stage-A Validation Report", f"", f"AUROC (z_max, pos vs neg): {au:.3f}", ""]
    for param, vals in [("rsasa_threshold", [0.15, 0.20, 0.25]),
                        ("patch_radius", [8.0, 10.0, 12.0]),
                        ("his_weight", [0.0, 0.05, 0.1])]:
        df = sweep_parameter(structs, cfg, param, vals)
        lines.append(f"- {param}: ranking stability (Spearman) = {ranking_stability(df):.3f}")
    rho = engine_agreement(structs, cfg)
    lines.append(f"- A-vs-B agreement (Spearman) = {rho:.3f}")
    Path("data/validation_report.md").write_text("\n".join(lines))
    print("\n".join(lines))

if __name__ == "__main__":
    main()
