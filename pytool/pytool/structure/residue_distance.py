import warnings

import click
import MDAnalysis as mda
from MDAnalysis.analysis import distances

warnings.filterwarnings("ignore")


def residue_distance(
    topology: str, residue1: int, residue2: int, dcd: str | None = None
) -> None:
    u = mda.Universe(topology) if dcd is None else mda.Universe(topology, dcd)

    res1 = u.select_atoms(f"name CA and resid {residue1}")
    res2 = u.select_atoms(f"name CA and resid {residue2}")

    if dcd is None:
        pass
    else:
        for _ts in u.trajectory:
            _, _, _dist = distances.dist(res1, res2)


@click.command()
@click.argument("topology", type=click.Path(exists=True), nargs=1)
@click.argument("residue1", type=click.INT, nargs=1)
@click.argument("residue2", type=click.INT, nargs=1)
@click.option("--dcd", default=None, type=click.Path(exists=True), nargs=1)
def residue_distance_to_command(
    topology: str, residue1: int, residue2: int, dcd: str | None
) -> None:
    residue_distance(topology, residue1, residue2, dcd)
