import argparse
from dit.model.scope import Scope
import os

def add_scope(directory: str):
    scope = Scope()
    assert os.path.isdir(directory), f"{directory} is not a directory"
    scope.add_directory(directory)

def to_command():
    parser = argparse.ArgumentParser()

    parser.add_argument("directory", type=str, help="新たに加える監視ディレクトリ")

    args = parser.parse_args()
    add_scope(args.directory)


