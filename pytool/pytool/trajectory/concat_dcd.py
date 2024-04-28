from glob import glob
import os
from typing import Optional
from argparse import ArgumentParser
from natsort import natsorted
from loguru import logger
import MDAnalysis as mda
import math

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
            logger.info(f"Found {len(dcd_files)} DCD files in {dcd_path}")
        elif os.path.isfile(dcd_path):
            result.append(dcd_path)
            logger.info(f"Found {dcd_path}")
        else:
            logger.warning(f"Cannot find dir/file {dcd_path}")
    return result

def concat_dcd(topology: str,dcds: list[str], output_name: str, dt_ratios: Optional[list[int]] = None):
    """Concatenate DCD files

    Args:
        dcds (list[str]): A list of DCD files
        output_name (str): Output file name
    """
    logger.info(f"Concatenating {len(dcds)} DCD files to {output_name}")

    u = mda.Universe(topology)

    assert u.atoms is not None, "No atoms found in the first DCD file"
    with mda.Writer(output_name, u.atoms.n_atoms) as W:
        if dt_ratios is None:
            u.load_new(dcds)
            for _ in u.trajectory:
                W.write(u)
        else:
            for dcd_path, dt_ratio in zip(dcds, dt_ratios):
                u.load_new(dcd_path)
                lcm = math.lcm(*dt_ratios)
                for i, _ in enumerate(u.trajectory):
                    if i % (lcm // dt_ratio) == 0:
                        W.write(u)


def main():
    parser = ArgumentParser()

    parser.add_argument("topology", type=str, help="Topology file")
    parser.add_argument("dcds", nargs="+", help="DCD files or dirs")
    parser.add_argument("-o", "--output", default="concat.dcd", help="Output file name. Default: concat.dcd")
    parser.add_argument("-r" ,"--dt-ratios", nargs="+", type=int, default=None,help="DT ratios for each DCD file. Default: None")

    args = parser.parse_args()
    logger.debug(args)

    dcd_paths = _glob_dcd(args.dcds)

    if len(dcd_paths) == 0:
        logger.error("No DCD files found")

    concat_dcd(args.topology, dcd_paths, args.output, args.dt_ratios)


if __name__ == "__main__":
    main()

