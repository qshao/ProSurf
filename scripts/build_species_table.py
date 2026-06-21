"""Build data/species_meltome_tms.csv from the Meltome Atlas workbooks.

One row per (protein, organism): uniprot, organism, ogt, organism_class,
tm (median across that organism's primary datasets), n_datasets.
Human is included here too; the scoring step reuses existing human scores.
"""
import csv
from pathlib import Path

from prosurf.io.meltome import species_datasets, organism_tms

ROOT = Path(__file__).resolve().parent.parent
S1 = ROOT / "data" / "meltome_s1.xlsx"
S4 = ROOT / "data" / "meltome_s4.xlsx"
OUT = ROOT / "data" / "species_meltome_tms.csv"


def main():
    spec = species_datasets(S1)
    print(f"Found {len(spec)} organisms in S1")
    rows = []
    for organism, info in sorted(spec.items(), key=lambda kv: kv[1]["ogt"]):
        tms = organism_tms(S4, info["dataset_ids"])
        print(f"  {organism:38s} OGT={info['ogt']:>4}  "
              f"datasets={len(info['dataset_ids'])}  proteins={len(tms)}")
        for uid, tm in tms.items():
            rows.append({
                "uniprot": uid,
                "organism": organism,
                "ogt": info["ogt"],
                "organism_class": info["class"],
                "tm": round(tm, 3),
                "n_datasets": len(info["dataset_ids"]),
            })
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=["uniprot", "organism", "ogt",
                           "organism_class", "tm", "n_datasets"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nWrote {len(rows)} rows to {OUT}")


if __name__ == "__main__":
    main()
