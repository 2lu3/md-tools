"""Analyze one equilibration run."""

from collections.abc import Callable
from dataclasses import dataclass
from typing import TYPE_CHECKING

import click
import matplotlib.pyplot as plt
from loguru import logger

from .common.plot import (
    plot_box_sizes,
    plot_pressure,
    plot_temperature,
    plot_total_energy,
)
from .common.reader import read_log

if TYPE_CHECKING:
    import pandas as pd


@dataclass
class EquilPlotParams:
    """Parameters used to plot an equilibration condition."""

    log_paths: list[str]
    condition_name: str


def analyze_equilibration(
    equil_list: list[EquilPlotParams],
    save_path_generator: Callable[[str], str],
    window_size: int = 10,
    *,
    popup: bool = False,
) -> None:
    """Analyze equilibration logs and save summary plots.

    Args:
        equil_list: Equilibration conditions to analyze.
        save_path_generator: Function mapping a feature name to an output path.
        window_size: Moving average window size.
        popup: Whether to show plots in a popup window.
    """
    df_list: list[pd.DataFrame] = [read_log(equil.log_paths) for equil in equil_list]

    condition_names: list[str] = [equil.condition_name for equil in equil_list]

    try:
        logger.info("Plotting Total Energy")
        plot_total_energy(
            df_list,
            condition_names,
            save_path_generator("total_energy"),
            window_size,
        )
    except KeyError as e:
        logger.warning(f"Could not plot Total Energy: {e}")
        logger.debug(e)

    try:
        logger.info("Plotting Temperature")
        plot_temperature(
            df_list,
            condition_names,
            save_path_generator("temperature"),
            window_size,
        )
    except KeyError as e:
        logger.warning(f"Could not plot Temperature: {e}")
        logger.debug(e)

    try:
        logger.info("Plotting Box Sizes")
        plot_box_sizes(df_list, condition_names, save_path_generator("box_sizes"))
    except KeyError as e:
        logger.warning(f"Could not plot Box Sizes: {e}")
        logger.debug(e)

    try:
        logger.info("Plotting Pressure")
        plot_pressure(
            df_list,
            condition_names,
            save_path_generator("pressure"),
            window_size,
        )
    except KeyError as e:
        logger.warning(f"Could not plot Pressure: {e}")
        logger.debug(e)

    if popup:
        plt.show()


@click.command()
@click.argument("log_path", type=click.Path(exists=True))
@click.option("--savename", type=str, default="equil", help="filename of the plot")
@click.option(
    "--window-size",
    type=int,
    default=10,
    help="window size for moving average",
)
@click.option("--popup", is_flag=True, help="show plots in popup window")
def command(
    log_path: str,
    savename: str,
    window_size: int,
    *,
    popup: bool,
) -> None:
    """Run equilibration analysis from the command line.

    Args:
        log_path: Log path passed by Click.
        savename: Output filename stem.
        window_size: Moving average window size.
        popup: Whether to show plots in a popup window.
    """
    analyze_equilibration(log_path, savename, window_size, popup=popup)
