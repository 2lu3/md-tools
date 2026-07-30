from functools import cache

import MDAnalysis as mda


def save_residue_selection(residue_selection: str, out_path: str) -> None:
    with open(out_path, "w") as f:
        f.write(residue_selection)


@cache
def load_residue_selection(residue_selection_path: str) -> str:
    with open(residue_selection_path) as f:
        return f.read()


def illustrate_selection(u: mda.Universe, selection: str, out_path: str) -> None:
    atoms = u.select_atoms(f"protein and {selection}")

    for atom in u.atoms:  # type: ignore
        atom.bfactor = 0

    for atom in atoms:
        atom.bfactor = 1

    u.atoms.write(out_path)  # type: ignore
