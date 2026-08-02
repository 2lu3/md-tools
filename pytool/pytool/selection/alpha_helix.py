"""Alpha-helix residue selection helpers."""

import warnings

import MDAnalysis as mda
from MDAnalysis.analysis.dssp import DSSP, translate

warnings.filterwarnings("ignore")

def _persistent_ss(u: mda.Universe, threshold: float) -> str:
    dssp = DSSP(u)
    long_run = dssp.run()
    persistent_residues = translate(
        long_run.results.dssp_ndarray.mean(axis=0) > threshold
    )
    return "".join(persistent_residues)



def _ss(u: mda.Universe) -> str:
    dssp = DSSP(u)
    long_run = dssp.run()
    return long_run.results.dssp[0]



def _merge_ss(ss_list: list[str], method: str) -> str:
    if not ss_list:
        msg = "ss_list must not be empty"
        raise ValueError(msg)

    # all ss should be same length
    ss_length_list = [len(ss) for ss in ss_list]
    if len(set(ss_length_list)) != 1:
        msg = "All secondary structures should be same length"
        raise ValueError(msg)
    ss_length = ss_length_list[0]

    merged_ss = ""
    if method == "all":
        for i in range(ss_length):
            if all(ss[i] == "H" for ss in ss_list):
                merged_ss += "H"
            else:
                merged_ss += "-"
    elif method == "any":
        for i in range(ss_length):
            if any(ss[i] == "H" for ss in ss_list):
                merged_ss += "H"
            else:
                merged_ss += "-"
    else:
        raise NotImplementedError

    return merged_ss


def select_alpha_helix(
    u_list: list[mda.Universe],
    method: str,
    threshold: float = 0.8,
) -> str:
    """Select residues that are persistently assigned as alpha helices.

    Args:
        u_list: Universes to inspect.
        method: Merge method for secondary-structure calls.
        threshold: DSSP occupancy threshold for persistent assignment.

    Returns:
        MDAnalysis selection string for alpha-helix residues.
    """
    if not u_list:
        msg = "u_list must not be empty"
        raise ValueError(msg)
    if not (0.0 <= threshold <= 1.0):
        msg = "threshold must be in [0.0, 1.0]"
        raise ValueError(msg)

    ss_list = []
    for u in u_list:
        ss = _persistent_ss(u, threshold)
        ss_list.append(ss)

    merged_ss = _merge_ss(ss_list, method=method)

    alpha_helix_residues = []
    u = u_list[0]
    # DSSP assigns states to protein residues, so align against protein-only order.
    protein_residues = u.select_atoms("protein").residues  # type: ignore[attr-defined]
    if len(protein_residues) != len(merged_ss):
        msg = "DSSP residue count does not match protein residues in the universe"
        raise ValueError(msg)

    for res, ss in zip(protein_residues, merged_ss, strict=False):
        if ss == "H":
            alpha_helix_residues.append(res.resid)

    if not alpha_helix_residues:
        msg = "No alpha-helix residues found"
        raise ValueError(msg)

    return f"(resid {' '.join(map(str, alpha_helix_residues))})"
