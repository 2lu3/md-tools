from collections.abc import Callable
from dataclasses import dataclass
from typing import TYPE_CHECKING

import click
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
    log_paths: list[str]
    condition_name: str

def analyze_equilibration(
    equil_list: list[EquilPlotParams],
    save_path_generator: Callable[[str], str],
    window_size=10,
    popup: bool = False,
) -> None:
    """analyze_equilibration.

    Args:
        equil_list (list[EquilPlotParams]): equil_list
        condition_names (list[str]): condition_names
        save_path_generator (Callable[[str], str]): df gen(feature_name) -> path
        window_size (int): window_size
        popup (bool): popup
    """
    df_list: list[pd.DataFrame] = []
    for equil in equil_list:
        df_list.append(read_log(equil.log_paths))

    condition_names: list[str] = [equil.condition_name for equil in equil_list]

    try:
        logger.info("Plotting Total Energy")
        plot_total_energy(df_list, condition_names, save_path_generator("total_energy"), window_size)
    except KeyError as e:
        logger.warning(f"Could not plot Total Energy: {e}")
        logger.debug(e)

    try:
        logger.info("Plotting Temperature")
        plot_temperature(df_list, condition_names, save_path_generator("temperature"), window_size)
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
        plot_pressure(df_list, condition_names, save_path_generator("pressure"), window_size)
    except KeyError as e:
        logger.warning(f"Could not plot Pressure: {e}")
        logger.debug(e)

    if popup:
        import matplotlib.pyplot as plt

        plt.show()


@click.command()
@click.argument("log_path", type=click.Path(exists=True))
@click.option("--savename", type=str, default="equil", help="filename of the plot")
@click.option("--window-size", type=int, default=10, help="window size for moving average")
@click.option("--popup", is_flag=True, help="show plots in popup window")
def command(log_path: str, savename: str, window_size: int, popup: bool) -> None:

    analyze_equilibration(log_path, savename, window_size, popup)
