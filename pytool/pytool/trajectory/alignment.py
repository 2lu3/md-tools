"""Align structures and trajectories."""

import warnings

import MDAnalysis as mda
from MDAnalysis.analysis.align import AlignTraj, alignto

warnings.filterwarnings("ignore")


def align(
    mobile: mda.Universe,
    ref: mda.Universe,
    output_path: str,
    select: str = "protein and name CA",
) -> None:
    """Align a mobile structure to a reference structure.

    Args:
        mobile: Mobile universe.
        ref: Reference universe.
        output_path: Output structure path.
        select: Atom selection used for alignment.
    """
    alignto(mobile, ref, select=select, weights="mass")

    mobile.atoms.write(output_path)


def align_trajectory(
    mobile: mda.Universe,
    ref: mda.Universe,
    output_path: str,
    select: str = "protein and name CA",
    *,
    verbose: bool = True,
) -> None:
    """Align a mobile trajectory to a reference structure.

    Args:
        mobile: Mobile universe.
        ref: Reference universe.
        output_path: Output trajectory path.
        select: Atom selection used for alignment.
        verbose: Whether to print MDAnalysis progress output.
    """
    AlignTraj(mobile, ref, select=select, filename=output_path, verbose=verbose).run()
