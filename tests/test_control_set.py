# tests/test_control_set.py
import yaml
from pathlib import Path

def test_control_set_has_both_classes():
    data = yaml.safe_load(Path("configs/control_set.yaml").read_text())
    assert len(data["positives"]) >= 5
    assert len(data["negatives"]) >= 5
    # UniProt-like accessions
    for k in ("positives", "negatives"):
        for acc in data[k]:
            assert isinstance(acc, str) and 6 <= len(acc) <= 10
