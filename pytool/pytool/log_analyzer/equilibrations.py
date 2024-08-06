from argparse import ArgumentParser
from glob import glob
from itertools import cycle
import os
from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from natsort import natsorted
import numpy as np
from dataclasses import dataclass

from .common.log_glob import glob_log_files
from .common.reader import read_column_by_name, read_column_names

linestyles = ["-", "--", "-.", ":"]
linestyle_cycler = cycle(linestyles)


def analyze_single_column(
    project_dirs: list[str],
    column_name: str,
    display_name: str,
    use_moving_average: bool,
    window_size: int,
    figsize: tuple[int, int],
) -> Figure:
    @dataclass
    class Data:
        project_name: str
        time: dict[str, list[int]]
        value: dict[str, list[float]]

    data_list: list[Data] = []
    for project_dir in project_dirs:
        data = Data(
            project_name=os.path.basename(os.path.normpath(project_dir)),
            time={},
            value={},
        )
        for log_file in glob_log_files(project_dir):
            if column_name not in read_column_names([log_file]):
                continue
            log_name = os.path.basename(os.path.normpath(log_file))
            time = list(map(int, read_column_by_name([log_file], "TIME")))
            value = list(map(float, read_column_by_name([log_file], column_name)))

            data.time[log_name] = time
            data.value[log_name] = value
        data_list.append(data)

    file_names = set(
        [file_name for data in data_list for file_name in data.time.keys()]
    )
    file_names = natsorted(file_names)

    fig = plt.figure(figsize=figsize, constrained_layout=True)

    for i, file_name in enumerate(file_names):
        ax = fig.add_subplot(1, len(file_names), i + 1)
        cmap = cm.get_cmap("tab20c")
        ax.set_prop_cycle(color=[cmap(i) for i in np.linspace(0, 1, len(data_list) // len(linestyles) + 1)] * len(linestyles))

        for data in data_list:
            if file_name in data.time.keys():
                if use_moving_average:
                    moving_average_x = np.convolve(
                        data.time[file_name],
                        np.ones(window_size) / window_size,
                        mode="valid",
                    )
                    moving_average_y = np.convolve(
                        data.value[file_name],
                        np.ones(window_size) / window_size,
                        mode="valid",
                    )
                    ax.plot(
                        moving_average_x,
                        moving_average_y,
                        label=f"{data.project_name}",
                        style=next(linestyle_cycler),
                    )
                else:
                    ax.plot(
                        data.time[file_name],
                        data.value[file_name],
                        label=data.project_name,
                        style=next(linestyle_cycler),
                    )
        ax.legend()
        ax.set_title(file_name)
        ax.set_xlabel("Time (ps)")
        ax.set_ylabel(column_name)

    fig.suptitle(display_name)
    return fig


def analyze_box_sizes(
    project_dirs: list[str],
    use_moving_average: bool = False,
    window_size: int = 10,
    figsize: tuple[int, int] = (12, 6),
):
    @dataclass
    class Data:
        project_name: str
        time: dict[str, list[int]]
        x: dict[str, list[float]]
        y: dict[str, list[float]]
        z: dict[str, list[float]]

    data_list: list[Data] = []
    for project_dir in project_dirs:
        data = Data(
            project_name=os.path.basename(os.path.normpath(project_dir)),
            time={},
            x={},
            y={},
            z={},
        )

        for log_file in glob_log_files(project_dir):
            if "BOXX" not in read_column_names([log_file]):
                continue

            log_name = os.path.basename(os.path.normpath(log_file))
            time = list(map(int, read_column_by_name([log_file], "TIME")))
            x = list(map(float, read_column_by_name([log_file], "BOXX")))
            y = list(map(float, read_column_by_name([log_file], "BOXY")))
            z = list(map(float, read_column_by_name([log_file], "BOXZ")))

            data.time[log_name] = time
            data.x[log_name] = x
            data.y[log_name] = y
            data.z[log_name] = z

        data_list.append(data)

    file_names = set(
        [file_name for data in data_list for file_name in data.time.keys()]
    )
    file_names = natsorted(file_names)

    fig = plt.figure(figsize=figsize, constrained_layout=True)

    index = 1
    for dim in "XYZ":
        for file_name in file_names:
            ax = fig.add_subplot(3, len(file_names), index)
            index += 1
            cmap = cm.get_cmap("tab20c")
            ax.set_prop_cycle(color=[cmap(i) for i in np.linspace(0, 1, len(data_list) // len(linestyles) + 1)] * len(linestyles))
            ax.set_prop_cycle(line=linestyles)

            ax.set_title(f"{file_name} {dim}")
            ax.set_xlabel("Time (ps)")
            ax.set_ylabel("Box Size (nm)")

            for data in data_list:
                if file_name in data.time.keys():
                    if dim == "X":
                        y = data.x
                    elif dim == "Y":
                        y = data.y
                    elif dim == "Z":
                        y = data.z
                    else:
                        raise ValueError(f"Invalid dimension {dim}")

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
                            style=next(linestyle_cycler),
                        )
                    else:
                        ax.plot(
                            data.time[file_name],
                            y[file_name],
                            label=data.project_name,
                            style=next(linestyle_cycler),
                        )

            ax.legend()

    fig.legend()
    fig.suptitle("Box Size")

    return fig


def analyze_equilibrations(
    project_dirs: list[str],
    use_moving_average: bool = False,
    window_size: int = 10,
    figsize: tuple[int, int] = (12, 6),
) -> list[tuple[str, Figure]]:
    """analyze_equilibrations.

    以下のグラフを作成する
    * 全エネルギー
    * ポテンシャルエネルギー
    * 運動エネルギー
    * 温度
    * 格子サイズ
    * 圧力


    Args:
        project_dirs (list[str]): project_dirs
        use_moving_average (bool): use_moving_average
        window_size (int): window_size
        figsize (tuple[int, int]): figsize

    Returns:
        Figure:
    """

    fig_total_energy = analyze_single_column(
        project_dirs,
        "TOTAL_ENE",
        "Total Energy (KJ/mol)",
        use_moving_average,
        window_size,
        figsize,
    )
    fig_potential_energy = analyze_single_column(
        project_dirs,
        "POTENTIAL_ENE",
        "Potential Energy (KJ/mol)",
        use_moving_average,
        window_size,
        figsize,
    )
    fig_kinetic_energy = analyze_single_column(
        project_dirs,
        "KINETIC_ENE",
        "Kinetic Energy (KJ/mol)",
        use_moving_average,
        window_size,
        figsize,
    )
    fig_temperature = analyze_single_column(
        project_dirs,
        "TEMPERATURE",
        "Temperature (K)",
        use_moving_average,
        window_size,
        figsize,
    )

    fig_box_sizes = analyze_box_sizes(
        project_dirs, use_moving_average, window_size, figsize
    )

    fig_pressure = analyze_single_column(
        project_dirs,
        "PRESSURE",
        "Pressure (bar)",
        use_moving_average,
        window_size,
        figsize,
    )

    return [
        ("total_energy", fig_total_energy),
        ("potential_energy", fig_potential_energy),
        ("kinetic_energy", fig_kinetic_energy),
        ("temperature", fig_temperature),
        ("box_sizes", fig_box_sizes),
        ("pressure", fig_pressure),
    ]


def command():
    parser = ArgumentParser()
    parser.add_argument("--project_dirs", nargs="+", type=str, default=None)
    parser.add_argument("--root-dir", type=str, default=None)
    parser.add_argument("--figsize", type=int, nargs=2, default=[12, 6])
    parser.add_argument("--out-name", type=str, default="equilibrations")
    parser.add_argument("--window-size", type=int, default=10)
    parser.add_argument("--use-moving-average", action="store_true")
    args = parser.parse_args()

    assert args.project_dirs is not None or args.root_dir is not None

    if args.project_dirs is None:
        dirs = glob(os.path.join(args.root_dir, "*/"))
    else:
        dirs = args.project_dirs

    figs = analyze_equilibrations(
        dirs,
        args.use_moving_average,
        args.window_size,
        args.figsize,
    )
    for name, fig in figs:
        fig.savefig(f"{args.out_name}_{name}.png")
