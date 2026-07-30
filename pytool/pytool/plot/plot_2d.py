
"""Plot selected columns from text data as 2D lines."""

from collections.abc import Sequence
from pathlib import Path

import click
import matplotlib.pyplot as plt


def plot_2d(  # noqa: PLR0913, PLR0917 - plotting API maps directly to CLI inputs.
    filename: str,
    x_index: int,
    y_indexes: Sequence[int],
    title: str,
    xlabel: str,
    ylabel: str,
    save_path: str | None = None,
) -> None:
    """Plot text-file columns as 2D lines.

    Args:
        filename: Input text file path.
        x_index: Column index to use for x values.
        y_indexes: Column indexes to plot as y values.
        title: Plot title.
        xlabel: X-axis label.
        ylabel: Y-axis label.
        save_path: Optional path to save the plot.
    """
    with Path(filename).open() as f:
        lines = f.readlines()
        x = []
        ys = [[] for _ in range(len(y_indexes))]
        for raw_line in lines:
            stripped_line = raw_line.strip()
            if stripped_line.startswith("#"):
                continue
            values = stripped_line.split()
            x.append(float(values[x_index]))
            for i, y_index in enumerate(y_indexes):
                ys[i].append(float(values[y_index]))
        for i, y in enumerate(ys):
            plt.plot(x, y, label=f"column {y_indexes[i]}")
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.legend()

        if save_path is not None:
            plt.savefig(save_path)
        else:
            plt.show()

@click.command()
@click.argument("filename")
@click.argument("x_index", type=int)
@click.argument("y_indexes", nargs=-1, type=int)
@click.option("--title", "-t", default="Plot", help="Title of the plot")
@click.option("--xlabel", "-x", default="x", help="Label of x axis")
@click.option("--ylabel", "-y", default="y", help="Label of y axis")
@click.option("--output_path", "-o", default=None, help="Path to save the plot")
def plot_2d_to_command(  # noqa: PLR0913, PLR0917 - click exposes each CLI input.
    filename: str,
    x_index: int,
    y_indexes: tuple[int, ...],
    title: str,
    xlabel: str,
    ylabel: str,
    output_path: str | None,
) -> None:
    """Run the plot-2d command."""
    plot_2d(filename, x_index, y_indexes, title, xlabel, ylabel, output_path)
