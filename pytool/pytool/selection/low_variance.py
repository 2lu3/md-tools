"""Low-variance residue selection helpers."""

import warnings

import MDAnalysis as mda
import numpy as np
from loguru import logger
from numpy.typing import NDArray

warnings.filterwarnings("ignore")

FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]
BoolArray = NDArray[np.bool_]


def load_position_resid(
    u_list: list[mda.Universe],
    selection: str | list[str],
) -> tuple[FloatArray, list[IntArray]]:
    """Load selected atom positions and residue IDs from universes.

    Args:
        u_list: Universes to read.
        selection: Atom selection string or per-universe selection strings.

    Returns:
        Tuple of positions and residue ID arrays.
    """
    logger.debug("Loading positions")

    selections = [selection] * len(u_list) if isinstance(selection, str) else selection

    positions: list[FloatArray] = []
    resids: list[IntArray] = []
    for u, sel in zip(u_list, selections, strict=False):
        for _ in u.trajectory:
            atoms = u.select_atoms(sel)
            positions.append(atoms.positions)
            resids.append(atoms.resids)

    positions = np.array(positions)
    return positions, resids


def kabsch(P: FloatArray, Q: FloatArray) -> tuple[FloatArray, FloatArray]:  # noqa: N803
    """Compute the optimal rotation and translation from P onto Q.

    Args:
        P: Moving point cloud with shape (N, 3).
        Q: Reference point cloud with shape (N, 3).

    Returns:
        Rotation matrix and translation vector.
    """
    centroid_p = np.mean(P, axis=0)
    centroid_q = np.mean(Q, axis=0)
    p_centered = P - centroid_p
    q_centered = Q - centroid_q
    covariance = np.dot(p_centered.T, q_centered)
    u_matrix, _singular_values, vt_matrix = np.linalg.svd(covariance)
    v_matrix = vt_matrix.T
    determinant = np.linalg.det(np.dot(v_matrix, u_matrix.T))
    correction = np.eye(3)
    correction[2, 2] = determinant
    rotation = np.dot(v_matrix, np.dot(correction, u_matrix.T))
    translation = centroid_q - np.dot(centroid_p, rotation)
    return rotation, translation


def align_positions(positions: FloatArray) -> FloatArray:
    """Align all coordinate frames to the mean structure.

    Args:
        positions: Coordinates with shape (n_frames, n_residues, 3).

    Returns:
        Aligned coordinates with the same shape.
    """
    n_frames = positions.shape[0]
    aligned_positions = np.empty_like(positions)
    ref = np.mean(positions, axis=0)

    for i in range(n_frames):
        frame_positions = positions[i]
        rotation, translation = kabsch(frame_positions, ref)
        aligned_positions[i] = np.dot(frame_positions, rotation) + translation
    return aligned_positions


def compute_rmsf(aligned_positions: FloatArray) -> FloatArray:
    """Compute per-residue RMSF values.

    Args:
        aligned_positions: Coordinates with shape (n_frames, n_residues, 3).

    Returns:
        RMSF values with shape (n_residues,).
    """
    mean_positions = np.mean(aligned_positions, axis=0)
    diffs = aligned_positions - mean_positions
    squared_diffs = np.sum(diffs**2, axis=2)
    return np.sqrt(np.mean(squared_diffs, axis=0))



def calc_low_variance(positions: FloatArray) -> tuple[list[float], list[BoolArray]]:
    """Calculate low-variance candidate masks by iterative RMSF pruning.

    Args:
        positions: Coordinates with shape (n_frames, n_residues, 3).

    Returns:
        Average RMSF values and candidate masks for each pruning step.
    """
    candidate_mask_list: list[BoolArray] = []
    avg_rmsf_list: list[float] = []

    candidate_mask = np.ones(positions.shape[1], dtype=bool)
    while candidate_mask.sum() > 0:
        positions_sel = positions[
            :, candidate_mask, :
        ]
        aligned_sel = align_positions(positions_sel)
        rmsf_sel = compute_rmsf(aligned_sel)
        avg_rmsf = np.mean(rmsf_sel)

        current_candidate_count = candidate_mask.sum()

        logger.info(
            f"Candidates: {current_candidate_count}, Average RMSF: {avg_rmsf:.3f} Å"
        )

        candidate_mask_list.append(candidate_mask.copy())
        avg_rmsf_list.append(avg_rmsf)

        idx_to_remove = np.argmax(rmsf_sel)
        candidate_indices = np.where(candidate_mask)[0]
        remove_index = candidate_indices[idx_to_remove]
        candidate_mask[remove_index] = False

    return avg_rmsf_list, candidate_mask_list


def calc_low_variance_residues(
    u_list: list[mda.Universe], selection: str | list[str] = "protein and name CA"
) -> list[tuple[float, list[str]]]:
    """Calculate low-variance residue IDs for each pruning step.

    Args:
        u_list: Universes to analyze.
        selection: Atom selection string or per-universe selection strings.

    Returns:
        Average RMSF values paired with remaining residue IDs.
    """
    positions, resids = load_position_resid(u_list, selection)

    if len({len(resid) for resid in resids}) != 1:
        msg = "All residues should have same length"
        raise ValueError(msg)

    avg_rmsfs, candidate_masks = calc_low_variance(positions)

    result: list[tuple[float, list[str]]] = []
    for candidate, rmsf in zip(candidate_masks, avg_rmsfs, strict=False):
        residues = list(map(str, resids[0][candidate]))
        result.append((rmsf, residues))

    return result


def select_low_variance_residues(
    u_list: list[mda.Universe],
    remain_res_num: int,
    selection: str | list[str] = "protein and name CA",
) -> str:
    """Select the requested count of lowest-variance residues.

    Args:
        u_list: Universes to analyze.
        remain_res_num: Number of residues to keep.
        selection: Atom selection string or per-universe selection strings.

    Returns:
        MDAnalysis selection string for low-variance residues.
    """
    result = calc_low_variance_residues(u_list, selection)

    residues = result[-remain_res_num][1]

    return "(resid " + " ".join(residues) + ")"
