import click
import subprocess

from dit.model.scope import Scope


def _git_pull(dry_run: bool):
    if dry_run:
        print("Dry run: git pull")
        return

    result = subprocess.run(["git", "pull"], capture_output=True, text=True)
    print(result)


def _dvc_pull(dry_run: bool):
    files_to_pull = []
    scope = Scope()
    scope.find_patterns(["*.dvc"])

    
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
    if dry_run:
        print("Dry run")
        return

    dvc_result = subprocess.run(["dvc", "pull"], capture_output=True, text=True)
    print(dvc_result)
    git_result = subprocess.run(["git", "pull"], capture_output=True, text=True)
    print(git_result)
