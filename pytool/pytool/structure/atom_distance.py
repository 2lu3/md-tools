"""Calculate distances between atom selections."""

import argparse
import warnings

import MDAnalysis as mda
from MDAnalysis.analysis import distances

warnings.filterwarnings("ignore")


def atom_distance(
    topology: str,
    sel1: str,
    sel2: str,
    dcd: str | list[str] | None = None,
) -> None:
    """Calculate distances between two atom selections.

    Args:
        topology: Topology file path.
        sel1: First atom selection.
        sel2: Second atom selection.
        dcd: Optional trajectory path or paths.
    """
    u = mda.Universe(topology) if dcd is None else mda.Universe(topology, dcd)

    selected1 = u.select_atoms(sel1)
    selected2 = u.select_atoms(sel2)

    if dcd is None:
        _, _, _dist = distances.dist(selected1, selected2)
    else:
        for _ts in u.trajectory:
            _, _, _dist = distances.dist(selected1, selected2)


def atom_distance_to_command() -> None:
    """Run the atom-distance command."""
    parser = argparse.ArgumentParser(description="Calculate distance between two atoms")
    parser.add_argument("topology", help="Topology file")
    parser.add_argument("sel1", help="Atom Selection 1")
    parser.add_argument("sel2", help="Atom Selection 2")
    parser.add_argument("--dcd", help="DCD files", default=None, nargs="+")
    args = parser.parse_args()
    atom_distance(args.topology, args.sel1, args.sel2, args.dcd)
