import MDAnalysis as mda
from functools import lru_cache


def save_residue_selection(residue_selection: str, out_path: str):
    with open(out_path, "w") as f:
        f.write(residue_selection)


@lru_cache(maxsize=None)
def load_residue_selection(residue_selection_path: str) -> str:
    with open(residue_selection_path, "r") as f:
        return f.read()


def illustrate_selection(u: mda.Universe, selection: str, out_path: str):
    atoms = u.select_atoms(f"protein and {selection}")

    for atom in u.atoms:  # type: ignore
        atom.bfactor = 0

    for atom in atoms:
        atom.bfactor = 1

    u.atoms.write(out_path)  # type: ignore
