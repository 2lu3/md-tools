from typing import Optional
from matplotlib import pyplot as plt

from .reader import read_column_by_name


def apply_window_size(x: list[float], window_size: int):
    result: list[float] = []
    for i in range(len(x) - window_size):
        result.append(sum(x[i : i + window_size]) / window_size)
    return result


def plot_2dline(
    log_paths: list[str],
    output_path: str,
    x_name: str,
    y_name: str,
    title: str,
    x_label: str,
    y_label: str,
    window_size: Optional[int]= None,
    use_y_log: bool = False,
):
    x = read_column_by_name(log_paths, x_name)
    y = read_column_by_name(log_paths, y_name)

    fig: plt.Figure = plt.figure()  # type: ignore
    ax: plt.Axes = fig.add_subplot(1, 1, 1)

    if window_size is not None:
        ax.plot(x, y, label="Raw")
        ax.plot(x[window_size:], apply_window_size(y, window_size), label="Moving average")
        ax.legend()
    else:
        ax.plot(x, y)

    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    if use_y_log:
        ax.set_yscale("log")

    fig.savefig(output_path)


def plot_box_sizes(log_paths: list[str], output_path: str, project_name: str):
    x = read_column_by_name(log_paths, "TIME")
    y_x = read_column_by_name(log_paths, "BOXX")
    y_y = read_column_by_name(log_paths, "BOXY")
    y_z = read_column_by_name(log_paths, "BOXZ")
    data_list = [
        ("X", y_x),
        ("Y", y_y),
        ("Z", y_z),
    ]

    fig: plt.Figure = plt.figure()  # type: ignore
    for i, data in enumerate(data_list):
        ax: plt.Axes = fig.add_subplot(3, 1, i + 1)
        ax.plot(x, data[1])
        ax.set_xlabel("Time (ps)")
        ax.set_ylabel(f"Box size (Å)")
        ax.set_title(f"{data[0]} Box size of {project_name}")

    fig.savefig(output_path)


def plot_pressure(log_paths: list[str], output_path: str, project_name: str, window_size: Optional[int]):
    plot_2dline(
        log_paths,
        output_path,
        "TIME",
        "PRESSURE",
        f"Pressure of {project_name}",
        "Time (ps)",
        "Pressure (bar)",
        window_size,
    )


def plot_temperature(log_paths: list[str], output_path: str, project_name: str, window_size: Optional[int]):
    plot_2dline(
        log_paths,
        output_path,
        "TIME",
        "TEMPERATURE",
        f"Temperature of {project_name}",
        "Time (ps)",
        "Temperature (K)",
        window_size,
    )


def plot_total_energy(log_paths: list[str], output_path: str, project_name: str, window_size: Optional[int], x_name: str = "TIME"):
    plot_2dline(
        log_paths,
        output_path,
        x_name,
        "TOTAL_ENE",
        f"Total Energy of {project_name}",
        "Time (ps)",
        "Energy (kJ/mol)",
        window_size,
    )
    
def plot_potential_energy(log_paths: list[str], output_path: str, project_name: str, window_size: Optional[int], use_y_log: bool = False,x_name: str = "TIME"):
    plot_2dline(
        log_paths,
        output_path,
        x_name,
        "POTENTIAL_ENE",
        f"Potential Energy of {project_name}",
        "Time (ps)",
        "Energy (kJ/mol)",
        window_size,
        use_y_log=use_y_log,
    )







