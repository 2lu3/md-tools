"""Copy bundled project templates into an output directory."""

import shutil
from pathlib import Path

import click
from loguru import logger


def copy_template(output_dir: str) -> None:
    """テンプレートのすべてのファイルをコピーする.

    mode: min, eq, pr
    """
    output_path = Path(output_dir)
    source_dir = _template_dir()
    target_dir = output_path / "template"

    shutil.copytree(source_dir, target_dir, dirs_exist_ok=True)
    logger.info(f"copy template files to {target_dir}")

    if not (output_path / "builder.py").exists():
        shutil.copy(target_dir / "builder.py", output_path)
        logger.info("copy builder.py to root direcotry")
    else:
        logger.warning("builder.py did not copy to root direcotry")


@click.command()
@click.option("--output_dir", "-o", default=".", help="output directory")
def copy_template_to_command(output_dir: str) -> None:
    """Run the copy-template command-line interface."""
    copy_template(output_dir)


def _template_dir() -> Path:
    return Path(__file__).parent.parent / "config" / "template"
