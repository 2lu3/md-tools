import click
import subprocess
import os

from dit.model.dvc import DVC


def _git_pull(dry_run: bool):
    if dry_run:
        print("Dry run: git pull")
        return

    proc = subprocess.Popen(
        ["git", "pull"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    while True:
        assert proc.stdout is not None
        line = proc.stdout.readline()

        if line:
            print(line)
            continue

        if proc.poll() is not None:
            break


def _dvc_pull(dry_run: bool):
    files_to_pull = []
    dvc = DVC()
    files_to_pull = dvc.find_dvc_files_in_scope()

    if dry_run:
        print("Dry run: dvc pull")
        print("Files to be pulled:")
        for file in files_to_pull:
            print(file)
        return

    proc = subprocess.Popen(
        ["dvc", "pull", *files_to_pull], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )

    while True:
        assert proc.stdout is not None
        line = proc.stdout.readline()

        if line:
            print(line)
            continue

        if proc.poll() is not None:
            break


@click.command()
@click.option("-d", "--dry-run", is_flag=True)
def pull(dry_run: bool):
    files_to_pull = DVC().find_dvc_files_in_scope()
    with open("tmp_execute.sh", "w") as f:
        f.write("#!/bin/bash\n")

        f.write("git pull\n")

        f.write("dvc pull ")
        for file in files_to_pull:
            f.write(file + " ")
        f.write("\n")

    if dry_run:
        subprocess.run(["cat", "tmp_execute.sh"])
    else:
        subprocess.run(["bash", "tmp_execute.sh"])
    os.remove("tmp_execute.sh")
