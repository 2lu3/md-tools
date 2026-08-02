"""Build comparison figures for a selected log column."""

from dataclasses import dataclass, field
from itertools import cycle
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.figure import Figure
from natsort import natsorted

from .log_glob import glob_log_files
from .reader import read_log

linestyles = ["-", "--", "-.", ":"]

@dataclass
class Data:
   """Log data grouped by project and file."""

   project_name: str
   time: dict[str, list[int]] = field(default_factory=dict)
   value: dict[str, list[float]] = field(default_factory=dict)


def read_column(project_dir: str, column_name: str) -> Data:
    """Read a selected column from all logs in a project.

    Args:
        project_dir: Directory or log file path to read.
        column_name: Column name to extract.

    Returns:
        Project data keyed by log file name.
    """
    data: Data = Data(
        project_name=Path(project_dir).name,
        time={},
        value={},
    )
    for log_file in glob_log_files(project_dir):
        df = read_log([log_file])
        if column_name not in df.columns:
            continue
        log_name = Path(log_file).name
        time = list(map(int, df["TIME"]))
        value = list(map(float, df[column_name]))

        data.time[log_name] = time
        data.value[log_name] = value
    return data

def read_projects_log(project_dirs: list[str], column_name: str) -> list[Data]:
    """Read a selected column from multiple projects.

    Args:
        project_dirs: Project directories or log files to read.
        column_name: Column name to extract.

    Returns:
        Naturally sorted project data.
    """
    data_list: list[Data] = [
        read_column(project_dir, column_name) for project_dir in project_dirs
    ]
    return natsorted(data_list, key=lambda x: x.project_name)

def fig_by_column(  # noqa: PLR0913
        project_dirs: list[str],
        column_name: str,
        display_name: str,
        *,
        use_moving_average: bool,
        window_size: int,
        figsize: tuple[int, int],
        )-> Figure:
    """Create a figure comparing one log column across projects.

    Args:
        project_dirs: Project directories or log files to read.
        column_name: Column name to plot.
        display_name: Figure title.
        use_moving_average: Whether to plot a moving average.
        window_size: Moving average window size.
        figsize: Figure size in inches.

    Returns:
        Comparison figure.
    """
    data_list = read_projects_log(project_dirs, column_name)


    file_names = {
        file_name for data in data_list for file_name in data.time
    }
    file_names = natsorted(file_names)

    fig = plt.figure(figsize=figsize, constrained_layout=True)

    for i, file_name in enumerate(file_names):
        ax = fig.add_subplot(1, len(file_names), i + 1)
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

        for data in data_list:
            if file_name in data.time:
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
                        linestyle=next(linestyle_cycler),
                    )
                else:
                    ax.plot(
                        data.time[file_name],
                        data.value[file_name],
                        label=data.project_name,
                        linestyle=next(linestyle_cycler),
                    )
        ax.legend()
        ax.set_title(file_name)
        ax.set_xlabel("Time (ps)")
        ax.set_ylabel(column_name)

    fig.suptitle(display_name)
    return fig

