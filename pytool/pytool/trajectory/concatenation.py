"""Concatenate DCD trajectory files."""

import warnings
from pathlib import Path

import MDAnalysis as mda
from loguru import logger
from natsort import natsorted

warnings.filterwarnings("ignore")


def _glob_dcd(dcds: list[str]) -> list[str]:
    """Glob DCD files from a list of dirs and files.

    Args:
        dcds (list[str]): A list of dirs and files.

    Returns:
        list[str]: A list of DCD files.

    Examples:
        >>> glob_dcd(["/path/to/dir1", "/path/to/dir2", "/path/to/file.dcd"])
        [
            "/path/to/dir1/1.dcd",
            "/path/to/dir1/2.dcd",
            "/path/to/dir2/1.dcd",
            "/path/to/dir2/2.dcd",
            "/path/to/file.dcd",
        ]
    """
    result: list[str] = []
    for dcd_path in dcds:
        path = Path(dcd_path)
        if path.is_dir():
            dcd_files = natsorted(str(file_path) for file_path in path.rglob("*.dcd"))
            result.extend(natsorted(dcd_files))
            logger.debug(f"Found {len(dcd_files)} DCD files in {dcd_path}")
        elif path.is_file():
            result.append(dcd_path)
            logger.debug(f"Found {dcd_path}")
        else:
            logger.warning(f"Cannot find dir/file {dcd_path}")
    return result


def concat_dcd(
    topology: str,
    dcds: list[str],
    output_name: str,
) -> None:
    """Concatenate DCD files.

    Args:
        topology (str): Topology file
        dcds (list[str]): A list of DCD files or directories
        output_name (str): Output file name
    """
    dcd_paths = _glob_dcd(dcds)

    if len(dcd_paths) == 0:
        logger.error("No DCD files found")
        return

    logger.debug(f"Concatenating {len(dcd_paths)} DCD files to {output_name}")

    u = mda.Universe(topology)

    if u.atoms is None:
        msg = "No atoms found in the first DCD file."
        raise ValueError(msg)
    with mda.Writer(output_name, u.atoms.n_atoms) as writer:
        for dcd_path in dcd_paths:
            u.load_new(dcd_path)
            for _ in u.trajectory:
                writer.write(u.atoms)
