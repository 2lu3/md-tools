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
            for file in glob.glob(os.path.join(path, "**", "*"), recursive=True):
                files_to_consider.append(file)
    return files_to_consider

@click.command()
@click.argument("path", nargs=-1)
@click.option("-A", "--all", is_flag=True)
@click.option("-d", "--dry-run", is_flag=True)
def add(paths: tuple[str, ...], is_all: bool, dry_run: bool):
    files_to_consider = search_files(set(paths))
    if is_all:
        scope = Scope()
        files_to_consider.extend(
                search_files(scope.directories)
                )

    extension = Extension()
    files_by_git = []
    files_by_dvc = []

    for file in files_to_consider:
        if extension.is_match(file):
            files_by_dvc.append(file)
        else:
            files_by_git.append(file)

    if dry_run:
        print("Dry run")
        print("Files to be added to git:")
        for file in files_by_git:
            print(file)
        print("Files to be added to dvc:")
        for file in files_by_dvc:
            print(file)
        return

    git_result = subprocess.run(["git", "add", " ".join(files_by_git)], capture_output=True, text=True)
    print(git_result)
    dvc_result = subprocess.run(["dvc", "add", " ".join(files_by_dvc)], capture_output=True, text=True)
    print(dvc_result)

