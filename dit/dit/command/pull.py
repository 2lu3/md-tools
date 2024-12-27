import click
import subprocess

from dit.model.dvc import DVC


def _git_pull(dry_run: bool):
    if dry_run:
        print("Dry run: git pull")
        return

    result = subprocess.run(["git", "pull"], capture_output=True, text=True)
    print(result)


def _dvc_pull(dry_run: bool):
    files_to_pull = []
    dvc = DVC()
    files_to_pull = dvc.find_dvc_files()

    if dry_run:
        print("Dry run: dvc pull")
        print("Files to be pulled:")
        for file in files_to_pull:
            print(file)
        return

    result = subprocess.run(["dvc", "pull", *files_to_pull], capture_output=True, text=True)
    print(result)

@click.command()
@click.option("-d", "--dry-run", is_flag=True)
def pull(dry_run: bool):
    _dvc_pull(dry_run)
    _git_pull(dry_run)
