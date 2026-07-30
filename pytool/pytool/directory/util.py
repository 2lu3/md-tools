"""Utilities for safely copying files into directory layouts."""

import shutil
from pathlib import Path

from loguru import logger


def copy_file_safe(
    source_file: str | Path,
    output_dir: str | Path,
    subdir: str,
    dest_filename: str,
) -> None:
    """Copy file after checking source file exists, creating output directory.

    Args:
        source_file: Source file path.
        output_dir: Output directory path.
        subdir: Subdirectory name under the output directory.
        dest_filename: Destination file name.
    """
    source_path = Path(source_file)
    if not source_path.exists():
        msg = f"{source_file} not found"
        raise FileNotFoundError(msg)

    dest_dir = Path(output_dir) / subdir
    dest_dir.mkdir(parents=True, exist_ok=True)

    dest_file = dest_dir / dest_filename
    shutil.copy(source_path, dest_file)
    logger.debug(f"Copy {source_file} to {dest_file}")


