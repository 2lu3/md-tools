
"""Calculate molecular structure box sizes."""

import click
import MDAnalysis as mda


def get_box_size(u: mda.Universe | str) -> tuple[float, float, float]:
    """Get box size from a universe.

    Args:
        u (mda.Universe): A universe

    Returns:
        tuple[float, float, float]: Box size in x, y, z
    """

    def get_box_size_from_atoms(atoms: mda.AtomGroup) -> tuple[float, float, float]:
        """Get box size from atoms.

        Args:
            atoms (mda.AtomGroup): Atoms

        Returns:
            list[float]: Box size in x, y, z
        """
        max_x = max(atoms.positions[:, 0])
        max_y = max(atoms.positions[:, 1])
        max_z = max(atoms.positions[:, 2])
        min_x = min(atoms.positions[:, 0])
        min_y = min(atoms.positions[:, 1])
        min_z = min(atoms.positions[:, 2])

        return (max_x - min_x, max_y - min_y, max_z - min_z)

    if isinstance(u, str):
        u = mda.Universe(u)

    if u.atoms is None:
        msg = "No atoms found in the first DCD file."
        raise ValueError(msg)
    return get_box_size_from_atoms(u.atoms)


@click.command()
@click.argument("input_pdb", type=click.Path(exists=True))
def get_box_size_to_command(input_pdb: str) -> None:
    """Run the get-boxsize command."""
    u = mda.Universe(input_pdb)

    get_box_size(u)
