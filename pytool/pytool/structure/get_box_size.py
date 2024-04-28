from typing import Union
import MDAnalysis as mda
import click

def get_box_size(u: Union[mda.Universe, str]) -> tuple[float, float, float]:
    """Get box size from a universe

    Args:
        u (mda.Universe): A universe

    Returns:
        tuple[float, float, float]: Box size in x, y, z
    """
    def get_box_size_from_atoms(atoms: mda.AtomGroup):
        """Get box size from atoms

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

    assert u.atoms is not None, "No atoms found in the first DCD file"
    return get_box_size_from_atoms(u.atoms)

@click.command()
@click.argument("input_pdb", type=click.Path(exists=True))
def get_box_size_to_command(input_pdb: str):
    u = mda.Universe(input_pdb)

    boxsize = get_box_size(u)

    print('dimensions:',boxsize[0])
    print('atoms:',boxsize[1])
