import click
import os
from dit.model.extension import Extension
from dit.model.scope import Scope
import glob
import subprocess
from loguru import logger

def search_files(path_list: set[str]):
    files_to_consider = []
    for path in path_list:
        path = os.path.normpath(path)

        if not os.path.exists(path):
            logger.error(f"{path} does not exist")
            raise FileNotFoundError

        if os.path.isfile(path):
            files_to_consider.append(path)
        else:
            files_to_consider.extend(
                glob.glob(os.path.join(path, "**", "*"), recursive=True)
                )
    return files_to_consider

@click.command()
@click.argument("paths", nargs=-1)
@click.option("-A", "--all-files", is_flag=True)
@click.option("-d", "--dry-run", is_flag=True)
def add(paths: tuple[str, ...], all_files: bool, dry_run: bool):
    files_to_consider = search_files(set(paths))
    if all_files:
        scope = Scope()
        files_to_consider.extend( search_files(scope.directories))

    extension = Extension()
    files_by_git = []
    files_by_dvc = []

    for file in files_to_consider:
        if extension.is_match(file):
            files_by_dvc.append(file)
        else:
            files_by_git.append(file)

    with open("tmp_execute.sh", "w") as f:
        f.write("#!/bin/bash\n")
        if files_by_git:
            f.write(f"git add {' '.join(files_by_git)}\n")
        if files_by_dvc:
            f.write(f"dvc add {' '.join(files_by_dvc)}\n")

    if dry_run:
        subprocess.run(["cat", "tmp_execute.sh"])
    else:
        subprocess.run(["bash", "tmp_execute.sh"])
    os.remove("tmp_execute.sh")
