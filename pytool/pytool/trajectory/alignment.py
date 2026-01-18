from typing import Optional, Union
import MDAnalysis as mda
from MDAnalysis.analysis.align import AlignTraj, alignto

import warnings

warnings.filterwarnings("ignore")


def align(
    mobile: str,
    ref: str,
    output_path: str,
    select: str = "protein and name CA",
):
    alignto(mda.Universe(mobile), mda.Universe(ref), select=select, weights="mass")

    output = mda.Writer(output_path, mobile.atoms.n_atoms)  # type: ignore
    output.write(mobile.atoms)


def align_trajectory(
    mobile: mda.Universe,
    ref: mda.Universe,
    output_path: str,
    select: str = "protein and name CA",
    verbose: bool = True,
):
    AlignTraj(mobile, ref, select=select, filename=output_path, verbose=verbose).run()
