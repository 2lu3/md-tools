import MDAnalysis as mda
from MDAnalysis.analysis.dssp import DSSP, translate
import warnings
import japanize_matplotlib
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")

def _persistent_ss(u: mda.Universe, threshold):
    dssp = DSSP(u)
    long_run = dssp.run()
    persistent_residues = translate(
        long_run.results.dssp_ndarray.mean(axis=0) > threshold
    )
    persistent_ss = "".join(persistent_residues)

    return persistent_ss


def _ss(u: mda.Universe):
    dssp = DSSP(u)
    long_run = dssp.run()
    ss = long_run.results.dssp[0]

    return ss


def _merge_ss(ss_list: list[str], method: str):
    if not ss_list:
        raise ValueError("ss_list must not be empty")

    # all ss should be same length
    ss_length_list = [len(ss) for ss in ss_list]
    if len(set(ss_length_list)) != 1:
        raise ValueError("All secondary structures should be same length")
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


def select_alpha_helix(u_list: list[mda.Universe], method: str,threshold: float = 0.8):
    if not u_list:
        raise ValueError("u_list must not be empty")
    if not (0.0 <= threshold <= 1.0):
        raise ValueError("threshold must be in [0.0, 1.0]")

    ss_list = []
    for u in u_list:
        ss = _persistent_ss(u, threshold)
        ss_list.append(ss)

    merged_ss = _merge_ss(ss_list, method=method)

    alpha_helix_residues = []
    u = u_list[0]
    # DSSP assigns states to protein residues, so align against protein-only order.
    protein_residues = u.select_atoms("protein").residues  # type: ignore
    if len(protein_residues) != len(merged_ss):
        raise ValueError(
            "DSSP residue count does not match protein residues in the universe"
        )

    for res, ss in zip(protein_residues, merged_ss):
        if ss == "H":
            alpha_helix_residues.append(res.resid)

    if not alpha_helix_residues:
        raise ValueError("No alpha-helix residues found")

    selection = f"(resid {' '.join(map(str, alpha_helix_residues))})"
    return selection
