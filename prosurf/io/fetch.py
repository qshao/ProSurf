from pathlib import Path
import requests

_API_BASE = "https://alphafold.ebi.ac.uk/api/prediction"
_FILE_BASE = "https://alphafold.ebi.ac.uk/files"


def af2_url(uniprot_id: str) -> str:
    """Resolve the current AF2 PDB URL for a UniProt ID via the EBI API."""
    r = requests.get(f"{_API_BASE}/{uniprot_id}", timeout=30)
    r.raise_for_status()
    entries = r.json()
    if not entries:
        raise ValueError(f"No AF2 prediction found for {uniprot_id}")
    return entries[0]["pdbUrl"]


def fetch_af2(uniprot_id: str, out_dir: Path) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    # Check for any cached version regardless of model version number
    existing = list(out_dir.glob(f"AF-{uniprot_id}-F1-model_v*.pdb"))
    if existing:
        return existing[0]
    url = af2_url(uniprot_id)
    fname = url.split("/")[-1]
    dest = out_dir / fname
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    dest.write_bytes(resp.content)
    return dest
