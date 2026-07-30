"""Analyze one minimization run."""

import click
import matplotlib.pyplot as plt
from loguru import logger

from .common.log_glob import glob_log_files
from .common.plot import plot_potential_energy


def analyze_minimization(
    log_path: str,
    filename: str = "min",
    title: str = "Minimization",
    window_size: int = 10,
    *,
    popup: bool = False,
) -> None:
    """Analyze minimization logs and save the energy plot.

    Args:
        log_path: Directory or log file path to analyze.
        filename: Output filename stem.
        title: Plot title.
        window_size: Moving average window size.
        popup: Whether to show plots in a popup window.
    """
    log_files = glob_log_files(log_path)

    logger.info(f"Found {len(log_files)} log files")

    plot_potential_energy(log_files, f"energy_{filename}", title, window_size, "STEP")
    logger.info("Saved Energy")

    logger.info("Saved all plots")

    if popup:
        plt.show()


@click.command()
@click.argument("log_path", type=click.Path(exists=True))
@click.option("--filename", type=str, default="energy", help="filename of the plot")
@click.option("--title", type=str, default="Minimization", help="title of the plot")
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
    """Run minimization analysis from the command line.

    Args:
        log_path: Log path passed by Click.
        filename: Output filename stem.
        title: Plot title.
        window_size: Moving average window size.
        popup: Whether to show plots in a popup window.
    """
    analyze_minimization(log_path, filename, title, window_size, popup=popup)
