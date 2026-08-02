"""Selection persistence and visualization helpers."""

from functools import cache
from pathlib import Path

import MDAnalysis as mda


def save_residue_selection(residue_selection: str, out_path: str) -> None:
    """Save a residue selection string.

    Args:
        residue_selection: Selection text to write.
        out_path: Output file path.
    """
    with Path(out_path).open("w") as f:
        f.write(residue_selection)


@cache
def load_residue_selection(residue_selection_path: str) -> str:
    """Load a residue selection string.

    Args:
        residue_selection_path: Input file path.

    Returns:
        Loaded selection text.
    """
    with Path(residue_selection_path).open() as f:
        return f.read()


def illustrate_selection(u: mda.Universe, selection: str, out_path: str) -> None:
    """Write a structure with selected residue B-factors set.

    Args:
        u: Universe containing atoms to annotate.
        selection: Residue selection expression.
        out_path: Output structure path.
    """
    atoms = u.select_atoms(f"protein and {selection}")

    for atom in u.atoms:  # type: ignore[attr-defined]
        atom.bfactor = 0

    for atom in atoms:
        atom.bfactor = 1

    u.atoms.write(out_path)  # type: ignore[attr-defined]
