from MDAnalysis.analysis.rms import RMSF
import MDAnalysis as mda
from MDAnalysis.analysis.align import AlignTraj
import MDAnalysis.transformations as trans



def rmsf(u: mda.Universe, selection: str):
    selected= u.select_atoms(selection)
    not_selected = u.select_atoms(f"not ({selection})")

    transformation = [
            trans.unwrap(selected),
            trans.center_in_box(selected, wrap=True),
            trans.wrap(not_selected),
            ]

    u.trajectory.add_transformations(*transformation)

    _= AlignTraj(u, u,select=selection, in_memory=True).run()

    ref_coords = u.trajectory.timeseries(asel=selection).mean(axis=1)

    reference = mda.Merge(selected).load_new(ref_coords[:, None, :], order="afc")

    rmsf = RMSF(selected, verbose=True).run()

    return rmsf.results.rmsf
