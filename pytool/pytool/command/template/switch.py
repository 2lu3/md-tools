#!/usr/env/bin python3

import glob
from argparse import ArgumentParser
import subprocess
from loguru import logger

def glob_genesis_dirs():
    return glob.glob("genesis*/")

def main(model: str):
    dirs = glob_genesis_dirs()

    for d in dirs:
        if not d.endswith(model):
            logger.debug(f"Skip {d}")
            continue

        subprocess.run(["unlink", "genesis"])
        subprocess.run(["ln", "-s", d, "genesis"])

        logger.info(f"Switched to {d}")

if __name__ == "__main__":
    parser = ArgumentParser("GENESISのGPUとCPUの切り替えを行う")
    parser.add_argument("model", type=str, help="cpu / gpu")
    args = parser.parse_args()
    main(args.model)

