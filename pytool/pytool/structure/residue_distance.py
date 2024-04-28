from typing import Optional
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import click

import warnings
warnings.filterwarnings("ignore")


def residue_distance(
    topology: str, residue1: int, residue2: int, dcd: Optional[str] = None
):
    if dcd is None:
        u = mda.Universe(topology)
    else:
        u = mda.Universe(topology, dcd)

    res1 = u.select_atoms(f"name CA and resid {residue1}")
    res2 = u.select_atoms(f"name CA and resid {residue2}")

    if dcd is None:
        print(distances.dist(res1, res2)[2][0])
    else:
        for ts in u.trajectory:
            _, _, dist = distances.dist(res1, res2)
            print(ts.time, dist[0])


@click.command()
@click.argument("topology", type=click.Path(exists=True), nargs=1)
@click.argument("residue1", type=click.INT, nargs=1)
@click.argument("residue2", type=click.INT, nargs=1)
@click.option("--dcd", default=None, type=click.Path(exists=True), nargs=1)
def residue_distance_to_command(
    topology: str, residue1: int, residue2: int, dcd: Optional[str]
):
    residue_distance(topology, residue1, residue2, dcd)
