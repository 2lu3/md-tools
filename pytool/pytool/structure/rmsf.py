"""Calculate RMSF values for an atom selection."""

import MDAnalysis as mda
import MDAnalysis.transformations as trans
import numpy as np
from MDAnalysis.analysis.align import AlignTraj
from MDAnalysis.analysis.rms import RMSF
from numpy.typing import NDArray


def rmsf(u: mda.Universe, selection: str) -> NDArray[np.float64]:
    """Calculate RMSF values for a selection.

    Args:
        u: Universe with trajectory.
        selection: Atom selection for RMSF.

    Returns:
        RMSF values.
    """
    selected = u.select_atoms(selection)
    not_selected = u.select_atoms(f"not ({selection})")

    transformation = [
        trans.unwrap(selected),
        trans.center_in_box(selected, wrap=True),
        trans.wrap(not_selected),
    ]

    u.trajectory.add_transformations(*transformation)

    AlignTraj(u, u, select=selection, in_memory=True).run()

    ref_coords = u.trajectory.timeseries(asel=selection).mean(axis=1)

    mda.Merge(selected).load_new(ref_coords[:, None, :], order="afc")

    rmsf = RMSF(selected, verbose=True).run()

    return rmsf.results.rmsf
