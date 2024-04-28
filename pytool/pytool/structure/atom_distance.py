from typing import Optional
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import argparse

import warnings

warnings.filterwarnings("ignore")


def atom_distance(topology: str, sel1: str, sel2: str, dcd: Optional[str] = None):
    if dcd is None:
        u = mda.Universe(topology)
    else:
        u = mda.Universe(topology, dcd)

    selected1 = u.select_atoms(sel1)
    selected2 = u.select_atoms(sel2)

    if dcd is None:
        _, _, dist = distances.dist(selected1, selected2)
        print(dist[0])
    else:
        for ts in u.trajectory:
            _, _, dist = distances.dist(selected1, selected2)
            print(ts.time, dist[0])


def atom_distance_to_command():
    parser = argparse.ArgumentParser(description="Calculate distance between two atoms")
    parser.add_argument("topology", help="Topology file")
    parser.add_argument("sel1", help="Atom Selection 1")
    parser.add_argument("sel2", help="Atom Selection 2")
    parser.add_argument("--dcd", help="DCD files", default=None, nargs="+")
    args = parser.parse_args()
    atom_distance(args.topology, args.sel1, args.sel2, args.dcd)
