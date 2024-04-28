import click
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
    popup: bool = False,
):
    log_files = glob_log_files(log_path)

    logger.info(f"Found {len(log_files)} log files")

    try:
        plot_total_energy(log_files, f"total_energy_{filename}", title, window_size, "TIME")
        logger.info(f"Saved Total Energy")
    except ValueError as e:
        logger.warning(f"Could not plot Total Energy: {e}")
        logger.debug(e)

    try:
        plot_temperature(log_files, f"temperature_{filename}", title, window_size)
        logger.info(f"Saved Temperature")
    except ValueError as e:
        logger.warning(f"Could not plot Temperature: {e}")
        logger.debug(e)

    try:
        plot_box_sizes(log_files, f"box_sizes_{filename}", title)
        logger.info(f"Saved Box Sizes")
    except ValueError as e:
        logger.warning(f"Could not plot Box Sizes: {e}")
        logger.debug(e)

    try:
        plot_pressure(log_files, f"pressure_{filename}", title, window_size)
        logger.info(f"Saved Pressure")
    except ValueError as e:
        logger.warning(f"Could not plot Pressure: {e}")
        logger.debug(e)

    logger.info(f"Saved all plots")

    if popup:
        import matplotlib.pyplot as plt

        plt.show()



@click.command()
@click.argument("log_path", type=click.Path(exists=True))
@click.option("--filename", type=str, default="pr", help="filename of the plot")
@click.option("--title", type=str, default="Production", help="title of the plot")
@click.option(
    "--window-size", type=int, default=10, help="window size for moving average"
)
@click.option("--popup", is_flag=True, help="show plots in popup window")
def command(log_path: str, filename: str, title: str, window_size: int, popup: bool):
    analyze_production(log_path, filename, title, window_size, popup)
