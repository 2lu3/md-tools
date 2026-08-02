"""Plot common log analysis figures."""

import pandas as pd
from matplotlib import pyplot as plt


def plot_box_sizes(
    df_list: list[pd.DataFrame], condition_names: list[str], output_path: str
) -> None:
    """Plot box size values for multiple conditions.

    Args:
        df_list: Data frames containing box size columns.
        condition_names: Labels for each condition.
        output_path: Path where the figure is saved.
    """
    fig = plt.figure()

    for i, (label, column_name) in enumerate(
        [("X", "BOXX"), ("Y", "BOXY"), ("Z", "BOXZ")]
    ):
        ax = fig.add_subplot(3, 1, i + 1)
        ax.set_xlabel("Time (ps)")
        ax.set_ylabel("Box size (Å)")
        ax.set_title(f"Box size {label}")

        for condition_name, df in zip(condition_names, df_list, strict=False):
            x = df["TIME"]
            y = df[column_name]
            ax.plot(x, y, label=condition_name)

        ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
        fig.tight_layout()

    fig.savefig(output_path)


def plot_pressure(
    df_list: list[pd.DataFrame],
    condition_names: list[str],
    output_path: str,
    window_size: int | None = None,
) -> None:
    """Plot pressure values for multiple conditions.

    Args:
        df_list: Data frames containing pressure values.
        condition_names: Labels for each condition.
        output_path: Path where the figure is saved.
        window_size: Optional moving average window size.
    """
    fig, ax = plt.subplots()
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Pressure (bar)")
    ax.set_title("Pressure")

    for condition_name, df in zip(condition_names, df_list, strict=False):
        time = (
            df["TIME"]
            if window_size is None
            else _apply_window_size(df["TIME"], window_size)
        )
        pressure = (
            df["PRESSURE"]
            if window_size is None
            else _apply_window_size(df["PRESSURE"], window_size)
        )
        ax.plot(time, pressure, label=condition_name)

    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    fig.tight_layout()
    fig.savefig(output_path)


def plot_temperature(
    df_list: list[pd.DataFrame],
    condition_names: list[str],
    output_path: str,
    window_size: int | None = None,
) -> None:
    """Plot temperature values for multiple conditions.

    Args:
        df_list: Data frames containing temperature values.
        condition_names: Labels for each condition.
        output_path: Path where the figure is saved.
        window_size: Optional moving average window size.
    """
    fig, ax = plt.subplots()
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Temperature (K)")
    ax.set_title("Temperature")

    for condition_name, df in zip(condition_names, df_list, strict=False):
        time = (
            df["TIME"]
            if window_size is None
            else _apply_window_size(df["TIME"], window_size)
        )
        temperature = (
            df["TEMPERATURE"]
            if window_size is None
            else _apply_window_size(df["TEMPERATURE"], window_size)
        )
        ax.plot(time, temperature, label=condition_name)

    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    fig.tight_layout()
    fig.savefig(output_path)


def plot_potential_energy(
    df_list: list[pd.DataFrame],
    condition_names: list[str],
    output_path: str,
    window_size: int | None = None,
) -> None:
    """Plot potential energy values for multiple conditions.

    Args:
        df_list: Data frames containing potential energy values.
        condition_names: Labels for each condition.
        output_path: Path where the figure is saved.
        window_size: Optional moving average window size.
    """
    fig, ax = plt.subplots()
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Potential Energy (kJ/mol)")
    ax.set_title("Potential Energy")

    for condition_name, df in zip(condition_names, df_list, strict=False):
        time = (
            df["TIME"]
            if window_size is None
            else _apply_window_size(df["TIME"], window_size)
        )
        potential_energy = (
            df["POTENTIAL_ENE"]
            if window_size is None
            else _apply_window_size(df["POTENTIAL_ENE"], window_size)
        )
        ax.plot(time, potential_energy, label=condition_name)

    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    fig.tight_layout()
    fig.savefig(output_path)


def plot_total_energy(
    df_list: list[pd.DataFrame],
    condition_names: list[str],
    output_path: str,
    window_size: int | None = None,
) -> None:
    """Plot total energy values for multiple conditions.

    Args:
        df_list: Data frames containing total energy values.
        condition_names: Labels for each condition.
        output_path: Path where the figure is saved.
        window_size: Optional moving average window size.
    """
    fig, ax = plt.subplots()
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Total Energy (kJ/mol)")
    ax.set_title("Total Energy")

    for condition_name, df in zip(condition_names, df_list, strict=False):
        total_energy = (
            df["TOTAL_ENE"]
            if window_size is None
            else _apply_window_size(df["TOTAL_ENE"], window_size)
        )
        time = (
            df["TIME"]
            if window_size is None
            else _apply_window_size(df["TIME"], window_size)
        )
        ax.plot(time, total_energy, label=condition_name)

    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    fig.tight_layout()
    fig.savefig(output_path)


def _apply_window_size(x: list[float], window_size: int) -> list[float]:
    return [
        sum(x[i : i + window_size]) / window_size
        for i in range(len(x) - window_size)
    ]
