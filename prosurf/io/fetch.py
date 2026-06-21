from pathlib import Path
import requests


def af2_url(uniprot_id: str) -> str:
    return f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"


def fetch_af2(uniprot_id: str, out_dir: Path) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / f"AF-{uniprot_id}-F1-model_v4.pdb"
    if dest.exists():
        return dest
    resp = requests.get(af2_url(uniprot_id), timeout=60)
    resp.raise_for_status()
    dest.write_bytes(resp.content)
    return dest
