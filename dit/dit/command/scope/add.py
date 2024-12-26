import argparse
from dit.model.scope import Scope
import os

def add_scope(args):
    scope = Scope()
    assert os.path.isdir(args.directory), f"{args.directory} is not a directory"
    scope.add_directory(args.directory)

def register_subparser(subparser):
    parser = subparser.add_parser("add", help="Scopeにディレクトリを追加する")
    parser.add_argument("directory", type=str, help="新たに加える監視ディレクトリ")

def handle(args: argparse.Namespace):
    if args.command != "add":
        return
    add_scope(args.directory)
