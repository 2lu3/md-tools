"""Run common DVC commands."""

import shutil
import subprocess
from collections.abc import Sequence

import click


def dvc_add(additional_globs: Sequence[str] | None = None) -> None:
    """Add common simulation artifacts to DVC.

    Args:
        additional_globs: Additional glob patterns to add.
    """
    if additional_globs is None:
        additional_globs = []
    dvc_executable = shutil.which("dvc")
    if dvc_executable is None:
        msg = "dvc executable was not found."
        raise FileNotFoundError(msg)
    subprocess.run(  # noqa: S603 - executable is resolved from PATH for this CLI.
        [
            dvc_executable,
            "add",
            "--glob",
            "**/*.dcd",
            "**/*.rst",
            "**/*.dvl",
            "**/*.npy",
            "**/*.pkl",
            "**/*.tar",
            *additional_globs,
        ],
        check=True,
    )

@click.command()
@click.option("--additional_globs", "-a", multiple=True, help="additional globs to add")
def dvc_add_to_command(additional_globs: tuple[str, ...]) -> None:
    """Run the dvc-add command."""
    dvc_add(additional_globs)
