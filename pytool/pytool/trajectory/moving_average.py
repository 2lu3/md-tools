import MDAnalysis as mda
import numpy as np
import argparse
from loguru import logger


def moving_average(pdb: str, dcd: str, window: int, out_dcd: str):
    logger.info(f"Calculating moving average with window size {window}")
    u = mda.Universe(pdb, dcd)

    # 全フレームの座標を取得
    positions = np.array([u.trajectory.ts.positions.copy() for _ in u.trajectory])

    # 移動平均を計算 (軸を指定して適用)
    # ここでは、軸0がフレーム、軸1が原子、軸2が座標(x, y, z)
    smoothed_positions = np.apply_along_axis(
        lambda m: np.convolve(m, np.ones(window) / window, mode="valid"),
        axis=0,
        arr=positions,
    )
    logger.debug(f"{positions.shape[0]} -> {smoothed_positions.shape[0]}")

    with mda.Writer(out_dcd, n_atoms=u.atoms.n_atoms) as writer:  # type: ignore
        for frame_positions in smoothed_positions:
            u.trajectory.ts.positions = frame_positions
            writer.write(u.atoms)


def moving_average_to_command():
    parser = argparse.ArgumentParser(
        description="Calculate a moving average of a DCD file"
    )

    parser.add_argument("pdb", help="PDB file")
    parser.add_argument("dcd", help="DCD file")
    parser.add_argument("out", help="Output DCD file")
    parser.add_argument("--window", type=int, help="Window size", default=10)

    args = parser.parse_args()

    moving_average(args.pdb, args.dcd, args.window, args.out)
