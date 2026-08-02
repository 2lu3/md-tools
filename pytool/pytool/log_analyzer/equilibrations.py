"""Analyze equilibration runs across multiple projects."""

from argparse import ArgumentParser
from dataclasses import dataclass
from itertools import cycle
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure
from matplotlib.pyplot import cm
from natsort import natsorted

from .common.fig_by_column import fig_by_column
from .common.log_glob import glob_log_files
from .common.reader import read_log

linestyles = ["-", "--", "-.", ":"]


def analyze_box_sizes(  # noqa: C901, PLR0912
    project_dirs: list[str],
    *,
    use_moving_average: bool = False,
    window_size: int = 10,
    figsize: tuple[int, int] = (12, 6),
) -> Figure:
    """Plot box size values across multiple projects.

    Args:
        project_dirs: Project directories or log files to analyze.
        use_moving_average: Whether to plot a moving average.
        window_size: Moving average window size.
        figsize: Figure size in inches.

    Returns:
        Box size comparison figure.
    """
    @dataclass
    class Data:
        """Box size data grouped by project and log file."""

        project_name: str
        time: dict[str, list[int]]
        x: dict[str, list[float]]
        y: dict[str, list[float]]
        z: dict[str, list[float]]

    data_list: list[Data] = []
    for project_dir in project_dirs:
        data = Data(
            project_name=Path(project_dir).name,
            time={},
            x={},
            y={},
            z={},
        )

        for log_file in glob_log_files(project_dir):
            df = read_log([log_file])
            if "BOXX" not in df.columns:
                continue

            log_name = Path(log_file).name
            time = list(map(int, df["TIME"]))
            x = list(map(float, df["BOXX"]))
            y = list(map(float, df["BOXY"]))
            z = list(map(float, df["BOXZ"]))

            data.time[log_name] = time
            data.x[log_name] = x
            data.y[log_name] = y
            data.z[log_name] = z

        data_list.append(data)

    data_list = sorted(data_list, key=lambda x: x.project_name)

    file_names = {
        file_name for data in data_list for file_name in data.time
    }
    file_names = natsorted(file_names)

    fig = plt.figure(figsize=figsize, constrained_layout=True)

    index = 1
    for dim in "XYZ":
        for file_name in file_names:
            ax = fig.add_subplot(3, len(file_names), index)
            index += 1
            cmap = cm.get_cmap("tab20c")
            ax.set_prop_cycle(
                color=[
                    cmap(color_index)
                    for color_index in np.linspace(
                        0,
                        1,
                        len(data_list) // len(linestyles) + 1,
                    )
                    for _ in range(len(linestyles))
                ]
            )
            linestyle_cycler = cycle(linestyles)

            ax.set_title(f"{file_name} {dim}")
            ax.set_xlabel("Time (ps)")
            ax.set_ylabel("Box Size (nm)")

            for data in data_list:
                if file_name in data.time:
                    if dim == "X":
                        y = data.x
                    elif dim == "Y":
                        y = data.y
                    elif dim == "Z":
                        y = data.z
                    else:
                        msg = f"Invalid dimension {dim}"
                        raise ValueError(msg)

                    if use_moving_average:
                        moving_average_x = np.convolve(
                            data.time[file_name],
                            np.ones(window_size) / window_size,
                            mode="valid",
                        )
                        moving_average_y = np.convolve(
                            data.y[file_name],
                            np.ones(window_size) / window_size,
                            mode="valid",
                        )
                        ax.plot(
                            moving_average_x,
                            moving_average_y,
                            label=data.project_name,
                            linestyle=next(linestyle_cycler),
                        )
                    else:
                        ax.plot(
                            data.time[file_name],
                            y[file_name],
                            label=data.project_name,
                            linestyle=next(linestyle_cycler),
                        )

            ax.legend()

    fig.legend()
    fig.suptitle("Box Size")

    return fig


def analyze_equilibrations(
    project_dirs: list[str],
    *,
    use_moving_average: bool = False,
    window_size: int = 10,
    figsize: tuple[int, int] = (12, 6),
) -> list[tuple[str, Figure]]:
    """Plot equilibration values across multiple projects.

    Creates total energy, potential energy, kinetic energy, temperature,
    box size, and pressure figures.

    Args:
        project_dirs: Project directories or log files to analyze.
        use_moving_average: Whether to plot a moving average.
        window_size: Moving average window size.
        figsize: Figure size in inches.

    Returns:
        Named equilibration comparison figures.
    """
    fig_total_energy = fig_by_column(
        project_dirs,
        "TOTAL_ENE",
        "Total Energy (KJ/mol)",
        use_moving_average=use_moving_average,
        window_size=window_size,
        figsize=figsize,
    )
    fig_potential_energy = fig_by_column(
        project_dirs,
        "POTENTIAL_ENE",
        "Potential Energy (KJ/mol)",
        use_moving_average=use_moving_average,
        window_size=window_size,
        figsize=figsize,
    )
    fig_kinetic_energy = fig_by_column(
        project_dirs,
        "KINETIC_ENE",
        "Kinetic Energy (KJ/mol)",
        use_moving_average=use_moving_average,
        window_size=window_size,
        figsize=figsize,
    )
    fig_temperature = fig_by_column(
        project_dirs,
        "TEMPERATURE",
        "Temperature (K)",
        use_moving_average=use_moving_average,
        window_size=window_size,
        figsize=figsize,
    )

    fig_box_sizes = analyze_box_sizes(
        project_dirs,
        use_moving_average=use_moving_average,
        window_size=window_size,
        figsize=figsize,
    )

    fig_pressure = fig_by_column(
        project_dirs,
        "PRESSURE",
        "Pressure (bar)",
        use_moving_average=use_moving_average,
        window_size=window_size,
        figsize=figsize,
    )

    return [
        ("total_energy", fig_total_energy),
        ("potential_energy", fig_potential_energy),
        ("kinetic_energy", fig_kinetic_energy),
        ("temperature", fig_temperature),
        ("box_sizes", fig_box_sizes),
        ("pressure", fig_pressure),
    ]


def command() -> None:
    """Run multi-project equilibration analysis from the command line."""
    parser = ArgumentParser()
    parser.add_argument("--project_dirs", nargs="+", type=str, default=None)
    parser.add_argument("--root-dir", type=str, default=None)
    parser.add_argument("--figsize", type=int, nargs=2, default=[12, 6])
    parser.add_argument("--out-name", type=str, default="equilibrations")
    parser.add_argument("--window-size", type=int, default=10)
    parser.add_argument("--use-moving-average", action="store_true")
    args = parser.parse_args()

    if args.project_dirs is None and args.root_dir is None:
        msg = "Either project_dirs or root_dir must be provided"
        raise ValueError(msg)

    if args.project_dirs is None:
        dirs = [str(project_dir) for project_dir in Path(args.root_dir).glob("*/")]
    else:
        dirs = args.project_dirs

    figs = analyze_equilibrations(
        dirs,
        use_moving_average=args.use_moving_average,
        window_size=args.window_size,
        figsize=args.figsize,
    )
    for name, fig in figs:
        fig.savefig(f"{args.out_name}_{name}.png")
