from argparse import ArgumentParser
from loguru import logger
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import os
import numpy as np

from .common.log_glob import glob_log_files
from .common.reader import read_column_by_name


def analyze_minimizations(
    project_dirs: list[str], window_size: int = 10, figsize: tuple[int, int] = (12, 6)
) -> Figure:
    """各種グラフを複数のプロジェクトを統合して描画する

    Args:
        project_dirs (list[str]): project_dirs
        window_size (int): window_size

    Returns:
        Figure:
    """
    names = []
    steps: list[list[int]] = []
    potential_energies: list[list[float]] = []
    for project_dir in project_dirs:
        log_files = glob_log_files(project_dir)

        logger.info(f"Found {len(log_files)} log files")

        step = read_column_by_name(log_files, "STEP")
        potential_energy = read_column_by_name(log_files, "POTENTIAL_ENE")

        step = list(map(int, step))
        potential_energy = list(map(float, potential_energy))

        names.append(os.path.basename(os.path.normpath(project_dir)))
        steps.append(step)
        potential_energies.append(potential_energy)

    fig: Figure = plt.figure(figsize=figsize)
    ax_all = fig.add_subplot(121)

    for name, step, potential_energy in zip(names, steps, potential_energies):
        ax_all.plot(step, potential_energy, label=name)

        moving_average_x = np.convolve(
            step, np.ones(window_size) / window_size, mode="valid"
        )
        moving_average_y = np.convolve(
            potential_energy, np.ones(window_size) / window_size, mode="valid"
        )
        ax_all.plot(
            moving_average_x, moving_average_y, label=f"{name} (Moving average)"
        )

    ax_all.set_title("Full period")
    ax_all.set_xlabel("STEP")
    ax_all.set_ylabel("Potential Energy (KJ/mol)")
    ax_all.legend()

    ax_half = fig.add_subplot(122)
    for name, step, potential_energy in zip(names, steps, potential_energies):
        half_index = len(step) // 2
        step = step[half_index:]
        potential_energy = potential_energy[half_index:]

        ax_half.plot(step, potential_energy, label=name)

        moving_average_x = np.convolve(
            step, np.ones(window_size) / window_size, mode="valid"
        )
        moving_average_y = np.convolve(
            potential_energy, np.ones(window_size) / window_size, mode="valid"
        )
        ax_half.plot(
            moving_average_x, moving_average_y, label=f"{name} (Moving average)"
        )

    ax_half.set_title("Second half")
    ax_half.set_xlabel("STEP")
    ax_half.set_ylabel("Potential Energy (KJ/mol))")
    ax_half.legend()

    fig.suptitle("Potential Energies")

    return fig


def command():
    parser = ArgumentParser()
    parser.add_argument("project_dirs", nargs="+", type=str)
    parser.add_argument("--figsize", type=int, nargs=2, default=[12, 6])
    parser.add_argument("--out", type=str, default="minimizations.png")
    parser.add_argument("--window_size", type=int, default=10)
    parser.add_argument("--popup", action="store_true")
    args = parser.parse_args()
    fig = analyze_minimizations(
        args.project_dirs, args.window_size, tuple([args.figsize[0], args.figsize[1]])
    )
    fig.savefig(args.out)

    if args.popup:
        plt.show()
