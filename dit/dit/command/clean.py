import click

from dit.model.dvc import DVC

# dic gc -w --not-in-remote

@click.command()
def clean():
    dvc = DVC()
    dvc.clean()
