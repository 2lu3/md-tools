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
    alignto(mobile, ref, select=select, weights="mass")

    mobile.atoms.write(output_path)

def align_trajectory(
    mobile: mda.Universe,
    ref: mda.Universe,
    output_path: str,
    select: str = "protein and name CA",
    verbose: bool = True,
) -> None:
    AlignTraj(mobile, ref, select=select, filename=output_path, verbose=verbose).run()
