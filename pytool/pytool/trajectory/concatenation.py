from glob import glob
import os
from natsort import natsorted
from loguru import logger
import MDAnalysis as mda

import warnings

warnings.filterwarnings("ignore")


def _glob_dcd(dcds: list[str]):
    """Glob DCD files from a list of dirs and files

    Args:
        dcds (list[str]): A list of dirs and files

    Returns:
        list[str]: A list of DCD files

    Examples:
        >>> glob_dcd(["/path/to/dir1", "/path/to/dir2", "/path/to/file.dcd"])
        ["/path/to/dir1/1.dcd", "/path/to/dir1/2.dcd", "/path/to/dir2/1.dcd", "/path/to/dir2/2.dcd", "/path/to/file.dcd"]
    """
    result: list[str] = []
    for dcd_path in dcds:
        if os.path.isdir(dcd_path):
            dcd_files = glob(os.path.join(dcd_path, "**", "*.dcd"), recursive=True)
            result.extend(natsorted(dcd_files))
            logger.debug(f"Found {len(dcd_files)} DCD files in {dcd_path}")
        elif os.path.isfile(dcd_path):
            result.append(dcd_path)
            logger.debug(f"Found {dcd_path}")
        else:
            logger.warning(f"Cannot find dir/file {dcd_path}")
    return result

def concat_dcd(
    topology: str,
    dcds: list[str],
    output_name: str,
):
    """Concatenate DCD files

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

    assert u.atoms is not None, "No atoms found in the first DCD file"
    with mda.Writer(output_name, u.atoms.n_atoms) as W:
        for dcd_path in dcd_paths:
            u.load_new(dcd_path)
            for _ in u.trajectory:
                W.write(u.atoms)
