#!/usr/bin/env python3

from pytool.config.render_template import render2file
from pytool.directory import (
    copy_toppar_files,
    copy_structure_files,
    clean_directory,
)
from pytool.structure import get_box_size
import argparse
import os

def create_project(project_dir: str, is_clean: bool, extra_keywords: dict = {}):
    project_name = os.path.basename(os.path.normpath(project_dir))
    clean_directory(project_dir, is_clean)
    copy_toppar_files("../01_data/c36", project_dir)
    copy_structure_files(
        output_dir=project_dir,
        pdb_file=f"../03_minimization/projects/{project_name}/pdb/input.pdb",
        psf_file=f"../03_minimization/projects/{project_name}/psf/input.psf",
        rst_file=f"../03_minimization/projects/{project_name}/out/min0.rst",
    )

    keywords: dict = {
        "nsteps": 100000,
        "mode": "eq",
        "index_num": 3,
        "node": 1,
        "proc": 8,
        "thread": 3,
        "proc_per_node": 8,
        "elapse": "01:00:00",
        "user_id": os.environ["FUGAKU_USER_ID"],
    }
    keywords.update(extra_keywords)

    for i in range(keywords["index_num"]):
        keywords["index"] = i
        render2file(f"{project_dir}/inp/{keywords['mode']}{i}.inp", f"{keywords['mode']}0.inp", keywords)
        render2file(f"{project_dir}/job{i}.sh", f"job.sh", keywords)
        os.chmod(f"{project_dir}/job{i}.sh", 0o755)

    render2file(f"{project_dir}/submit.sh", "submit.sh", keywords)
    os.chmod(f"{project_dir}/submit.sh", 0o755)

def create_benchmark():
    nodes = [1, 2, 4, 8]
    proc_per_node = 24
    thread = 2
    for node in nodes:
        keywords = {
            "node": node,
            "proc_per_node": proc_per_node,
            "thread": thread,
            "elapse": "00:10:00",
            "proc": node * proc_per_node,
        }
        create_project(f"benchmarks/project_{node}", True, keywords)


def main(is_clean: bool):
    projects = []
    for project in projects:
        create_project(f"projects/{project}", is_clean)
    create_benchmark()

    render2file("submit_all.sh", "submit_all.sh", {})
    os.chmod("submit_all.sh", 0o755)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--clean", action="store_true")

    args = parser.parse_args()

    main(args.clean)
