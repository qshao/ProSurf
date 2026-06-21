"""Score every protein in data/species_meltome_tms.csv across all organisms.

Resumable: appends to data/multispecies_scores.csv and skips proteins already
present there or in data/multispecies_failed.txt. Length filter < 150 aa.

Flags:
  --organisms "A,B"   restrict to these organisms (substring match)
  --limit N           cap proteins scored per organism (smoke testing)
  --min-len N         length floor (default 150)
"""
import argparse
import csv
import time
from pathlib import Path

import biotite.structure as struc

from prosurf.config import MetricConfig, PathsConfig
from prosurf.io.fetch import fetch_af2
from prosurf.io.parse import load_structure
from prosurf.metric.composition import composition_features
from prosurf.pipeline import analyze_structure

ROOT = Path(__file__).resolve().parent.parent
IN_CSV = ROOT / "data" / "species_meltome_tms.csv"
OUT_CSV = ROOT / "data" / "multispecies_scores.csv"
FAILED = ROOT / "data" / "multispecies_failed.txt"

FIELDS = ["uniprot", "organism", "ogt", "tm", "z_mean", "z_max", "z_frac",
          "n_patches", "charged_frac", "surface_charged_frac",
          "net_charge_per_res", "surface_charge_density"]


def load_done():
    done = set()
    if OUT_CSV.exists():
        with open(OUT_CSV) as f:
            for r in csv.DictReader(f):
                done.add((r["uniprot"], r["organism"]))
    return done


def load_failed():
    if FAILED.exists():
        return {line for line in FAILED.read_text().splitlines() if line}
    return set()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--organisms", default="")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--min-len", type=int, default=150)
    args = ap.parse_args()

    cfg = MetricConfig()
    paths = PathsConfig()

    with open(IN_CSV) as f:
        targets = list(csv.DictReader(f))
    if args.organisms:
        wanted = [s.strip().lower() for s in args.organisms.split(",")]
        targets = [t for t in targets
                   if any(w in t["organism"].lower() for w in wanted)]

    done = load_done()
    failed = load_failed()

    # Seed existing human scores from the full-proteome run (reuse, don't refetch).
    human_seed = ROOT / "data" / "all_9160_scores.csv"
    seeded = {}
    if human_seed.exists():
        with open(human_seed) as f:
            for r in csv.DictReader(f):
                seeded[r["uniprot"]] = r

    if not OUT_CSV.exists():
        with open(OUT_CSV, "w", newline="") as f:
            csv.DictWriter(f, fieldnames=FIELDS).writeheader()

    # Per-organism caps for --limit.  Initialise from already-done rows so
    # a resumed run with the same --limit produces 0 new proteins.
    per_org_count: dict = {}
    if args.limit:
        with open(IN_CSV) as f:
            all_targets = list(csv.DictReader(f))
        for t in all_targets:
            if (t["uniprot"], t["organism"]) in done:
                org = t["organism"]
                per_org_count[org] = per_org_count.get(org, 0) + 1
    t0 = time.time()
    n_new = 0

    for t in targets:
        uid, organism = t["uniprot"], t["organism"]
        if (uid, organism) in done or f"{uid}:{organism}" in failed:
            continue
        if args.limit:
            c = per_org_count.get(organism, 0)
            if c >= args.limit:
                continue

        # For human proteins already scored in the full-proteome run, seed from
        # that file to avoid re-downloading and re-scoring.
        if organism == "Homo sapiens" and uid in seeded:
            s = seeded[uid]
            # Try to compute composition features from the cached PDB file.
            comp = {"charged_frac": float("nan"), "surface_charged_frac": float("nan"),
                    "net_charge_per_res": float("nan"), "surface_charge_density": float("nan")}
            existing_pdbs = list(paths.data_dir.glob(f"AF-{uid}-F1-model_v*.pdb"))
            if existing_pdbs:
                try:
                    arr = load_structure(existing_pdbs[0])
                    comp = composition_features(arr, cfg)
                    comp = {k: float(v) for k, v in comp.items()}
                except Exception:
                    pass
            row = {
                "uniprot": uid,
                "organism": organism,
                "ogt": t["ogt"],
                "tm": t["tm"],
                "z_mean": s["z_mean"],
                "z_max": s["z_max"],
                "z_frac": s["z_frac"],
                "n_patches": s["n_patches"],
                **comp,
            }
            with open(OUT_CSV, "a", newline="") as f:
                csv.DictWriter(f, fieldnames=FIELDS).writerow(row)
            n_new += 1
            per_org_count[organism] = per_org_count.get(organism, 0) + 1
            done.add((uid, organism))
            if n_new % 25 == 0:
                rate = n_new / (time.time() - t0)
                print(f"  {n_new} seeded/scored  ({organism} {uid})  "
                      f"{rate:.1f}/s", flush=True)
            continue

        try:
            pdb = fetch_af2(uid, paths.data_dir)
            arr = load_structure(pdb)
            n_res = len(struc.get_residues(arr)[0])
            if n_res < args.min_len:
                with open(FAILED, "a") as ff:
                    ff.write(f"{uid}:{organism}\n")
                continue
            _, ps = analyze_structure(pdb, uid, cfg)
            comp = composition_features(arr, cfg)
            row = {
                "uniprot": uid, "organism": organism, "ogt": t["ogt"], "tm": t["tm"],
                "z_mean": ps.z_mean, "z_max": ps.z_max, "z_frac": ps.z_frac,
                "n_patches": ps.n_patches, **comp,
            }
            with open(OUT_CSV, "a", newline="") as f:
                csv.DictWriter(f, fieldnames=FIELDS).writerow(row)
            n_new += 1
            per_org_count[organism] = per_org_count.get(organism, 0) + 1
            done.add((uid, organism))
            if n_new % 25 == 0:
                rate = n_new / (time.time() - t0)
                print(f"  {n_new} scored  ({organism} {uid} "
                      f"z_mean={ps.z_mean:.3f})  {rate:.1f}/s", flush=True)
        except Exception as e:
            with open(FAILED, "a") as ff:
                ff.write(f"{uid}:{organism}\n")
            print(f"  FAIL {organism} {uid}: {e}", flush=True)

    print(f"Done. {n_new} new proteins scored -> {OUT_CSV}")


if __name__ == "__main__":
    main()
