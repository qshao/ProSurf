from pathlib import Path
import biotite.structure as struc
from biotite.structure.io.pdb import PDBFile


def load_structure(path: Path) -> struc.AtomArray:
    pdb = PDBFile.read(str(path))
    arr = pdb.get_structure(model=1)
    arr = arr[struc.filter_amino_acids(arr)]
    return arr
