"""Analyze one production run."""

import click
import matplotlib.pyplot as plt
from loguru import logger

from .common.log_glob import glob_log_files
from .common.plot import (
    plot_box_sizes,
    plot_pressure,
    plot_temperature,
    plot_total_energy,
)


def analyze_production(
    log_path: str,
    filename: str = "pr",
    title: str = "Production",
    window_size: int = 10,
    *,
    popup: bool = False,
) -> None:
    """Analyze production logs and save summary plots.

    Args:
        log_path: Directory or log file path to analyze.
        filename: Output filename stem.
        title: Plot title.
        window_size: Moving average window size.
        popup: Whether to show plots in a popup window.
    """
    log_files = glob_log_files(log_path)

    logger.info(f"Found {len(log_files)} log files")

    try:
        plot_total_energy(
            log_files,
            f"total_energy_{filename}",
            title,
            window_size,
            "TIME",
        )
        logger.info("Saved Total Energy")
    except ValueError as e:
        logger.warning(f"Could not plot Total Energy: {e}")
        logger.debug(e)

    try:
        plot_temperature(log_files, f"temperature_{filename}", title, window_size)
        logger.info("Saved Temperature")
    except ValueError as e:
        logger.warning(f"Could not plot Temperature: {e}")
        logger.debug(e)

    try:
        plot_box_sizes(log_files, f"box_sizes_{filename}", title)
        logger.info("Saved Box Sizes")
    except ValueError as e:
        logger.warning(f"Could not plot Box Sizes: {e}")
        logger.debug(e)

    try:
        plot_pressure(log_files, f"pressure_{filename}", title, window_size)
        logger.info("Saved Pressure")
    except ValueError as e:
        logger.warning(f"Could not plot Pressure: {e}")
        logger.debug(e)

    logger.info("Saved all plots")

    if popup:
        plt.show()



@click.command()
@click.argument("log_path", type=click.Path(exists=True))
@click.option("--filename", type=str, default="pr", help="filename of the plot")
@click.option("--title", type=str, default="Production", help="title of the plot")
@click.option(
    "--window-size", type=int, default=10, help="window size for moving average"
)
@click.option("--popup", is_flag=True, help="show plots in popup window")
def command(
    log_path: str,
    filename: str,
    title: str,
    window_size: int,
    *,
    popup: bool,
) -> None:
    """Run production analysis from the command line.

    Args:
        log_path: Log path passed by Click.
        filename: Output filename stem.
        title: Plot title.
        window_size: Moving average window size.
        popup: Whether to show plots in a popup window.
    """
    analyze_production(log_path, filename, title, window_size, popup=popup)
