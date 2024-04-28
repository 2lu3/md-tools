#!/usr/bin/env python3
import re
import glob
import os
from dataclasses import dataclass, field
from typing import Optional
from loguru import logger


class Result:
    def __init__(self, dir_path: str, name: str, node: int, proc: int, thread: int):
        self.dir_path: str = dir_path
        self.name: str = name
        self.node: int = node
        self.proc: int = proc
        self.thread: int = thread

    def validate(self):
        if len(self._gather_log_files()) == 0:
            return "No log files found"
        try:
            self.dynamics_time
            self.nsteps
        except RuntimeError as e:
            return e

        return "valid"

    @property
    def dynamics_time(self):
        for log_file in self._gather_log_files():
            with open(log_file) as f:
                for line in f.readlines():
                    if "dynamics" in line:
                        return float(line.split("=")[1])
        raise RuntimeError("Could not find dynamics time")

    @property
    def nsteps(self):
        for log_file in self._gather_log_files():
            with open(log_file) as f:
                for line in f.readlines():
                    if "nsteps" in line:
                        return int(line.strip().split("=")[-1])
        raise RuntimeError("Could not find nsteps")

    @property
    def core_efficiency(self):
        return self.nstep_per_sec / self.proc / self.thread

    @property
    def node_efficiency(self):
        return self.nstep_per_sec / self.node

    @property
    def nstep_per_sec(self):
        return self.nsteps / self.dynamics_time




    def _gather_log_files(self):
        files = glob.glob(os.path.join(self.dir_path, "out", "*.log*"))
        return files


def glob_benchmarks(projects_root_path: str):
    results = []
    pattern_re = re.compile(r".*?\/([a-z0-9]+)_(\d{1,3})_(\d{1,3})_(\d{1,3}).*?")

    for directory in glob.glob(os.path.join(projects_root_path, "*/")):
        match = pattern_re.match(directory)
        if match:
            results.append(
                Result(
                    directory,
                    match.group(1),
                    int(match.group(2)),
                    int(match.group(3)),
                    int(match.group(4)),
                )
            )
        else:
            print(f"Skipping {directory} because it does not match the pattern")

    return results


def analyze_benchmark(projects_root_path: str):
    logger.info(f"benchmark dir: {projects_root_path}/benchmarks")

    results: list[Result] = []
    for result in glob_benchmarks(projects_root_path):
        if result.validate() != "valid":
            logger.warning(f"Invalid result {result.dir_path}: {result.validate()}")
            continue
        results.append(result)

    for name in set([result.name for result in results]):
        grouped_results = [result for result in results if result.name == name]

        grouped_results.sort(key=lambda x: x.core_efficiency, reverse=True)

        # print in markdown table format with aligned for bash output
        print(f"## {name}")
        print(
            "| node | proc | thread | total time | nstep/sec | eff(node) | eff(core) |"
        )
        print(
            "|------|------|--------|------------|-----------|-----------|-----------|"
        )
        for result in grouped_results:
            print(
                f"| {result.node:4d} | {result.proc:4d} | {result.thread:6d} | {result.dynamics_time:10.2f} | {result.nstep_per_sec:9.2f} | {result.node_efficiency:9.2f} | {result.core_efficiency:9.2f} |"
                    )

if __name__ == "__main__":
    analyze_benchmark("benchmarks")
