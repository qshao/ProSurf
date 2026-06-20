from pathlib import Path
from prosurf.io.fetch import af2_url
from prosurf.io.parse import load_structure
import biotite.structure as struc


def test_af2_url():
    assert af2_url("P69905") == "https://alphafold.ebi.ac.uk/files/AF-P69905-F1-model_v4.pdb"


def test_load_structure_amino_acids_only(tmp_path):
    arr = load_structure(Path("tests/data/mini.pdb"))
    res_ids, res_names = struc.get_residues(arr)
    assert set(res_names) == {"LYS", "ASP", "ALA"}
    assert len(res_ids) == 3
