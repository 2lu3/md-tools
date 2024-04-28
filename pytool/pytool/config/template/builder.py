#!/usr/bin/env python3

from pytool.config.project_builder import ProjectBuilder
from pytool.directory import (
    copy_toppar_files,
    clean_directory,
)
from pytool.structure import get_box_size
import shutil
import argparse
import os


class Builder(ProjectBuilder):
    def project_params(self):
        params_hoge: dict = {
            "project_name": self.project_name,  # 固定
            "project_dir": "projects/" + self.project_name, # 固定
            "pdb_path": f"../03_minimization/{self.project_name}/pdb/input.pdb",
            "psf_path": f"../03_minimization/{self.project_name}/psf/input.psf",
            "rst_path": f"../03_minimization/{self.project_name}/out/min0.rst",
        }

        params = [
            params_hoge,
        ]

        for param in params:
            param.update(self._common_params())
            yield param

    def benchmark_params(self):
        resources = [
            {"node": 1, "proc": 16, "thread": 3},
            {"node": 2, "proc": 32, "thread": 3},
            {"node": 4, "proc": 64, "thread": 3},
            {"node": 8, "proc": 128, "thread": 3},
        ]

        for resource in resources:
            for param in self.project_params():
                bench_params = {
                    "project_dir": f"benchmarks/{param['project_name']}_{resource['node']}_{resource['proc']}_{resource['thread']}",
                    "node": resource["node"],
                    "proc": resource["proc"],
                    "thread": resource["thread"],
                    "proc_per_node": resource["proc"] // resource["node"],
                    "elapse": "00:10:00",
                    "indexes": [0],
                    "nsteps": 1000,
                }
                param.update(bench_params)
                yield param

    def create_project(self, param):
        clean_directory(param["project_dir"], self.is_clean_dirs)

        copy_toppar_files(param["project_dir"], "../01_data/c36")

        # pdb
        shutil.copy(
            param["pdb_path"],
            f"{param['project_dir']}/pdb/input.pdb",
        )

        # psf
        shutil.copy(
            param["psf_path"],
            f"{param['project_dir']}/psf/input.psf",
        )

        # rst file
        shutil.copy(
            param["rst_path"],
            f"{param['project_dir']}/rst/input.rst",
        )

        for index in param["indexes"]:
            param["index"] = index
            self.write_input(**param)
            self.write_job(**param)
            self.write_submit(**param)

    def _common_params(self):
        """同じproject名から生成するプロジェクトで共通するパラーメーターを返す"""
        return {
            "mode": "eq",
            "box_size": (0, 0, 0),
            # "box_size": get_box_size(f"../01_data/{self.project_name}/genesis/step3_input.pdb"),
            "node": 1,
            "proc": 16,
            "thread": 3,
            "proc_per_node": 16,
            "elapse": "01:00:00",
            "indexes": [0, 1, 2],
            "nsteps": 100000,
            "append_submit_all": False,
            "user_id": os.environ["FUGAKU_USER_ID"],
        }


def main(is_clean: bool):
    projects = []
    for project in projects:
        builder = Builder(project, is_benchmark=False, is_clean_dirs=is_clean)
        builder.build()

        builder = Builder(project, is_benchmark=True, is_clean_dirs=is_clean)
        builder.build()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--clean", action="store_true")

    args = parser.parse_args()

    main(args.clean)
