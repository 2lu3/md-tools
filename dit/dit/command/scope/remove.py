import argparse
from dit.model.scope import Scope
import os

def remove_scope(directory: str):
    scope = Scope()
    scope.remove_directory(directory)

def to_command():
    parser = argparse.ArgumentParser()

    parser.add_argument("directory", type=str, help="削除する監視ディレクトリ")

    args = parser.parse_args()
    remove_scope(args.directory)


