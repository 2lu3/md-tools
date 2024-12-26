import argparse
from dit.model.scope import Scope
import os

def add_scope(directory: str):
    scope = Scope()
    assert os.path.isdir(directory), f"{directory} is not a directory"
    scope.add_directory(directory)

def register_subparser(subparser):
    parser = subparser.add_parser("add-scope", help="Scopeにディレクトリを追加する")
    parser.add_argument("directory", type=str, help="新たに加える監視ディレクトリ")

def handle(args: argparse.Namespace):
    add_scope(args.directory)
