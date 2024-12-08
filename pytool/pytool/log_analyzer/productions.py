from argparse import ArgumentParser
from glob import glob
import os
from matplotlib.figure import Figure

from .common.fig_by_column import fig_by_column


def analyze_productions(
    project_dirs: list[str],
    use_moving_average: bool = False,
    window_size: int = 10,
    figsize: tuple[int, int] = (12, 6),
    ) -> list[tuple[str, Figure]]:

    fig_total_energy = fig_by_column(
            project_dirs,
            "TOTAL_ENE",
            "Total Energy (KJ/mol)",
            use_moving_average,
            window_size,
            figsize,
            )

    fig_potential_energy = fig_by_column(
            project_dirs,
            "POTENTIAL_ENE",
            "Potential Energy (KJ/mol)",
            use_moving_average,
            window_size,
            figsize,
            )

    fig_kinetic_energy = fig_by_column(
            project_dirs,
            "KINETIC_ENE",
            "Kinetic Energy (KJ/mol)",
            use_moving_average,
            window_size,
            figsize,
            )

    fig_temperature = fig_by_column(
            project_dirs,
            "TEMPERATURE",
            "Temperature (K)",
            use_moving_average,
            window_size,
            figsize,
            )

    fig_pressure = fig_by_column(
            project_dirs,
            "PRESSURE",
            "Pressure (bar)",
            use_moving_average,
            window_size,
            figsize,
            )

    return [
        ("total_energy", fig_total_energy),
        ("potential_energy", fig_potential_energy),
        ("kinetic_energy", fig_kinetic_energy),
        ("temperature", fig_temperature),
        ("pressure", fig_pressure),
    ]


def command():
    parser = ArgumentParser()
    parser.add_argument("--project_dirs", nargs="+", type=str, default=None)
    parser.add_argument("--root-dir", type=str, default=None)
    parser.add_argument("--figsize", type=int, nargs=2, default=[12, 6])
    parser.add_argument("--out-name", type=str, default="productions")
    parser.add_argument("--window-size", type=int, default=10)
    parser.add_argument("--use-moving-average", action="store_true")

    args = parser.parse_args()

    assert args.project_dirs is not None or args.root_dir is not None, "Either project_dirs or root_dir must be provided"

    if args.project_dirs is None:
        dirs = glob(os.path.join(args.root_dir, "*/"))
    else:
        dirs = [os.path.abspath(project_dir) for project_dir in args.project_dirs]

    figs = analyze_productions(
        dirs,
        args.use_moving_average,
        args.window_size,
        args.figsize,
    )

    for name, fig in figs:
        fig.savefig(f"{args.out_name}_{name}.png")

