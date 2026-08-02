"""Analyze production runs across multiple projects."""

from argparse import ArgumentParser
from pathlib import Path

from matplotlib.figure import Figure

from .common.fig_by_column import fig_by_column


def analyze_productions(
    project_dirs: list[str],
    *,
    use_moving_average: bool = False,
    window_size: int = 10,
    figsize: tuple[int, int] = (12, 6),
) -> list[tuple[str, Figure]]:
    """Plot production values across multiple projects.

    Args:
        project_dirs: Project directories or log files to analyze.
        use_moving_average: Whether to plot a moving average.
        window_size: Moving average window size.
        figsize: Figure size in inches.

    Returns:
        Named production comparison figures.
    """
    fig_total_energy = fig_by_column(
        project_dirs,
        "TOTAL_ENE",
        "Total Energy (KJ/mol)",
        use_moving_average=use_moving_average,
        window_size=window_size,
        figsize=figsize,
    )

    fig_potential_energy = fig_by_column(
        project_dirs,
        "POTENTIAL_ENE",
        "Potential Energy (KJ/mol)",
        use_moving_average=use_moving_average,
        window_size=window_size,
        figsize=figsize,
    )

    fig_kinetic_energy = fig_by_column(
        project_dirs,
        "KINETIC_ENE",
        "Kinetic Energy (KJ/mol)",
        use_moving_average=use_moving_average,
        window_size=window_size,
        figsize=figsize,
    )

    fig_temperature = fig_by_column(
        project_dirs,
        "TEMPERATURE",
        "Temperature (K)",
        use_moving_average=use_moving_average,
        window_size=window_size,
        figsize=figsize,
    )

    fig_pressure = fig_by_column(
        project_dirs,
        "PRESSURE",
        "Pressure (bar)",
        use_moving_average=use_moving_average,
        window_size=window_size,
        figsize=figsize,
    )

    return [
        ("total_energy", fig_total_energy),
        ("potential_energy", fig_potential_energy),
        ("kinetic_energy", fig_kinetic_energy),
        ("temperature", fig_temperature),
        ("pressure", fig_pressure),
    ]


def command() -> None:
    """Run multi-project production analysis from the command line."""
    parser = ArgumentParser()
    parser.add_argument("--project_dirs", nargs="+", type=str, default=None)
    parser.add_argument("--root-dir", type=str, default=None)
    parser.add_argument("--figsize", type=int, nargs=2, default=[12, 6])
    parser.add_argument("--out-name", type=str, default="productions")
    parser.add_argument("--window-size", type=int, default=10)
    parser.add_argument("--use-moving-average", action="store_true")

    args = parser.parse_args()

    if args.project_dirs is None and args.root_dir is None:
        msg = "Either project_dirs or root_dir must be provided"
        raise ValueError(msg)

    if args.project_dirs is None:
        dirs = [str(project_dir) for project_dir in Path(args.root_dir).glob("*/")]
    else:
        dirs = [str(Path(project_dir).resolve()) for project_dir in args.project_dirs]

    figs = analyze_productions(
        dirs,
        use_moving_average=args.use_moving_average,
        window_size=args.window_size,
        figsize=args.figsize,
    )

    for name, fig in figs:
        fig.savefig(f"{args.out_name}_{name}.png")

