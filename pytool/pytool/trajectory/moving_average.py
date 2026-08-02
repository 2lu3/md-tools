"""Smooth DCD trajectories with a moving average."""

import argparse

import MDAnalysis as mda
import numpy as np
from loguru import logger


def moving_average(pdb: str, dcd: str, window: int, out_dcd: str) -> None:
    """Write a moving-average-smoothed DCD trajectory.

    Args:
        pdb: PDB topology file path.
        dcd: Input DCD trajectory path.
        window: Moving-average window size.
        out_dcd: Output DCD trajectory path.
    """
    logger.info(f"Calculating moving average with window size {window}")
    u = mda.Universe(pdb, dcd)

    positions = np.array([u.trajectory.ts.positions.copy() for _ in u.trajectory])

    smoothed_positions = np.apply_along_axis(
        lambda m: np.convolve(m, np.ones(window) / window, mode="valid"),
        axis=0,
        arr=positions,
    )
    logger.debug(f"{positions.shape[0]} -> {smoothed_positions.shape[0]}")

    with mda.Writer(out_dcd, n_atoms=u.atoms.n_atoms) as writer:
        for frame_positions in smoothed_positions:
            u.trajectory.ts.positions = frame_positions
            writer.write(u.atoms)


def moving_average_to_command() -> None:
    """Run the dcd-moving-average command."""
    parser = argparse.ArgumentParser(
        description="Calculate a moving average of a DCD file"
    )

    parser.add_argument("pdb", help="PDB file")
    parser.add_argument("dcd", help="DCD file")
    parser.add_argument("out", help="Output DCD file")
    parser.add_argument("--window", type=int, help="Window size", default=10)

    args = parser.parse_args()

    moving_average(args.pdb, args.dcd, args.window, args.out)
