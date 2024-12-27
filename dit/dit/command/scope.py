from dit.model.scope import Scope
import os
from natsort import natsorted
from loguru import logger
import click

@click.group()
def scope():
    pass

@scope.command()
@click.argument("directory")
def add(directory):
    scope = Scope()
    assert os.path.isdir(directory), f"{directory} is not a directory"
    scope.add(directory)
    logger.info(f"add {directory} to scope")

@scope.command()
def list():
    scope = Scope()
    directories = scope.directories
    for directory in natsorted(directories):
        print(directory)

@scope.command()
@click.argument("directory")
def remove(directory):
    scope = Scope()
    scope.remove(directory)
    logger.info(f"remove {directory} from scope")

