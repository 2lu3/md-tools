"""Increment numbered directory prefixes."""

import shutil
import sys
from pathlib import Path

import click
from loguru import logger
from natsort import natsorted


def _setup_logger(*, verbose: bool) -> None:
    if verbose:
        logger.remove()
        logger.add(
            sys.stderr,
            format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {function} | {message}",
        )
    else:
        logger.remove()


def _glob_numbered_dirs() -> list[str]:
    """Return numbered directories in natural sort order.

    i.e.
    .
    ├── 00_hello
    ├── 01_hoge
    ├── 02_fuga
    └── 03_piyo.
    """
    dirs = [path for path in Path().glob("[0-9][0-9]_*") if path.is_dir()]
    normalized_dirs = [str(path) for path in dirs]
    return natsorted(normalized_dirs)


def _increment_directory_index(dir_name: str, add_num: int) -> str:
    """Increment a directory number by the requested amount.

    i.e.
    03_piyo -> 04_piyo.
    """
    number = int(dir_name[:2])
    new_number = number + add_num
    new_dir = dir_name.replace(f"{number:02}", f"{new_number:02}")
    shutil.move(dir_name, new_dir)

    return new_dir


def increment_directory_index(
    start_dir: str,
    add_num: int = 1,
    verbose: bool = False,  # noqa: FBT001, FBT002
) -> None:
    """Increment numbered directories from the selected start directory.

    i.e.
    start_dir: 02_fuga
    00_hello -> 00_hello
    01_hoge -> 01_hoge
    02_fuga -> 03_fuga
    03_piyo -> 04_piyo.
    """
    _setup_logger(verbose=verbose)

    numbered_dirs = _glob_numbered_dirs()
    logger.info(f"Found {len(numbered_dirs)} directories. {numbered_dirs}")

    # get index of start_dir
    start_idx = numbered_dirs.index(str(Path(start_dir)))

    # increment index
    for dir_name in numbered_dirs[start_idx:]:
        new_dir = _increment_directory_index(dir_name, add_num)
        logger.info(f"{dir_name} -> {new_dir}")


@click.command()
@click.argument("start_dir", type=click.Path(exists=True))
@click.option("--add_num", "-a", default=1, help="add number")
@click.option("--verbose", "-v", is_flag=True, help="Verbose mode")
def increment_directory_index_to_command(
    start_dir: str,
    add_num: int,
    verbose: bool,  # noqa: FBT001
) -> None:
    """Run the increment-dir command-line interface."""
    increment_directory_index(start_dir, add_num, verbose)
