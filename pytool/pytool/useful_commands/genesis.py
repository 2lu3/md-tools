"""Install GENESIS helper scripts."""

import shutil
import subprocess
from pathlib import Path

import click


def _bash_executable() -> str:
    bash = shutil.which("bash")
    if bash is None:
        msg = "bash executable was not found."
        raise FileNotFoundError(msg)
    return bash


@click.command()
def install_genesis_to_command() -> None:
    """Run the GENESIS installer script."""
    script_path = Path(__file__).parent / "install_genesis.sh"

    subprocess.run(  # noqa: S603 - script path is bundled with this package.
        [_bash_executable(), str(script_path)],
        check=True,
    )


@click.command()
def install_genesis_requirements_to_command() -> None:
    """Run the GENESIS requirements installer script."""
    script_path = Path(__file__).parent / "install_genesis_requirements.sh"

    subprocess.run(  # noqa: S603 - script path is bundled with this package.
        [_bash_executable(), str(script_path)],
        check=True,
    )
