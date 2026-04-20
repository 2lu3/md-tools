from typing import Optional
import pandas as pd
from matplotlib import pyplot as plt

from .reader import read_column_by_name, read_log


def plot_box_sizes(
    df_list: list[pd.DataFrame], condition_names: list[str], output_path: str
):
    fig = plt.figure()

    for i, label, column_name in enumerate(
        [("X", "BOXX"), ("Y", "BOXY"), ("Z", "BOXZ")]
    ):
        ax = fig.add_subplot(3, 1, i + 1)
        ax.set_xlabel("Time (ps)")
        ax.set_ylabel(f"Box size (Å)")
        ax.set_title(f"Box size {label}")

        for df in df_list:
            x = df["TIME"]
            y = df[column_name]

            ax.plot(x, y, label=condition_names[i])

        ax.legend()

    fig.savefig(output_path)


def plot_pressure(
    df_list: list[pd.DataFrame],
    condition_names: list[str],
    output_path: str,
    window_size: Optional[int] = None,
):
    fig, ax = plt.subplots()
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Pressure (bar)")
    ax.set_title("Pressure")

    for condition_name in condition_names:
        for df in df_list:
            pressure = df["PRESSURE"] if window_size is None else _apply_window_size(df["PRESSURE"], window_size)
            ax.plot(df["TIME"], pressure, label=condition_name)

    ax.legend()
    fig.savefig(output_path)


def plot_temperature(
    df_list: list[pd.DataFrame],
    condition_names: list[str],
    output_path: str,
    window_size: Optional[int] = None,
):
    fig, ax = plt.subplots()
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Temperature (K)")
    ax.set_title("Temperature")

    for condition_name in condition_names:
        for df in df_list:
            temperature = df["TEMPERATURE"] if window_size is None else _apply_window_size(df["TEMPERATURE"], window_size)
            ax.plot(df["TIME"], temperature, label=condition_name)

    ax.legend()
    fig.savefig(output_path)

def plot_potential_energy(
    df_list: list[pd.DataFrame],
    condition_names: list[str],
    output_path: str,
    window_size: Optional[int] = None,
):

    fig, ax = plt.subplots()
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Potential Energy (kJ/mol)")
    ax.set_title("Potential Energy")

    for condition_name in condition_names:
        for df in df_list:
            potential_energy = df["POTENTIAL_ENE"] if window_size is None else _apply_window_size(df["POTENTIAL_ENE"], window_size)
            ax.plot(df["TIME"], potential_energy, label=condition_name)

    ax.legend()
    fig.savefig(output_path)


def _apply_window_size(x: list[float], window_size: int):
    result: list[float] = []
    for i in range(len(x) - window_size):
        result.append(sum(x[i : i + window_size]) / window_size)
    return result