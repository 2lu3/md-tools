import click
from loguru import logger

from .common.plot import plot_potential_energy
from .common.log_glob import glob_log_files


def analyze_minimization(log_path: str, filename: str="min",title: str="Minimization",window_size: int=10, popup: bool=False):
    log_files = glob_log_files(log_path)

    logger.info(f"Found {len(log_files)} log files")

    plot_potential_energy(log_files, f"energy_{filename}", title, window_size, "STEP")
    logger.info(f"Saved Energy")

    logger.info(f"Saved all plots")

    if popup:
        import matplotlib.pyplot as plt

        plt.show()


@click.command()
@click.argument("log_path", type=click.Path(exists=True))
@click.option("--filename", type=str, default="energy", help="filename of the plot")
@click.option("--title", type=str, default="Minimization", help="title of the plot")
@click.option(
    "--window-size", type=int, default=10, help="window size for moving average"
)
@click.option("--popup", is_flag=True, help="show plots in popup window")
def command(log_path: str, filename: str,title: str,window_size: int, popup: bool):
    analyze_minimization(log_path, filename,title,window_size, popup)
