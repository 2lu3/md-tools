import click
import os
from loguru import logger

from .common.plot import plot_box_sizes, plot_pressure, plot_temperature, plot_total_energy
from .common.log_glob import glob_log_files


def analyze_equilibration(
    log_path: str,
    savename: str = "equil",
    window_size=10,
    popup: bool = False,
):
    log_files = glob_log_files(log_path)

    logger.info(f"Found {len(log_files)} log files")

    for log_file in log_files:
        logger.info(f"Processing {log_file}")
        basename = os.path.basename(log_file).split(".")[0]

        try:
            plot_total_energy([log_file], f"total_energy_{savename}_{basename}", basename, window_size, "TIME")
            logger.info(f"Saved Total Energy")
        except ValueError as e:
            logger.warning(f"Could not plot Total Energy for {log_file}")
            logger.debug(e)

        try:
            plot_temperature([log_file], f"temperature_{savename}_{basename}", basename, window_size)
            logger.info(f"Saved Temperature")
        except ValueError as e:
            logger.warning(f"Could not plot Temperature for {log_file}")
            logger.debug(e)

        try:
            plot_box_sizes([log_file], f"box_sizes_{savename}_{basename}", basename)
            logger.info(f"Saved Box Sizes")
        except ValueError as e:
            logger.warning(f"Could not plot Box Sizes for {log_file}")
            logger.debug(e)

        try:
            plot_pressure([log_file], f"pressure_{savename}_{basename}", basename, window_size)
            logger.info(f"Saved Pressure")
        except ValueError as e:
            logger.warning(f"Could not plot Pressure for {log_file}")
            logger.debug(e)

    logger.info(f"Saved all plots")

    if popup:
        import matplotlib.pyplot as plt

        plt.show()


@click.command()
@click.argument("log_path", type=click.Path(exists=True))
@click.option("--savename", type=str, default="equil", help="filename of the plot")
@click.option("--window-size", type=int, default=10, help="window size for moving average")
@click.option("--popup", is_flag=True, help="show plots in popup window")
def command(log_path: str, savename: str, window_size: int, popup: bool):

    analyze_equilibration(log_path, savename, window_size, popup)
