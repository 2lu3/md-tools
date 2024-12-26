import argparse
from dit.model.scope import Scope


def remove_scope(directory: str):
    scope = Scope()
    scope.remove_directory(directory)


def register_subparser(subparser):
    parser = subparser.add_parser(
        "remove", help="Scopeからディレクトリを削除する"
    )
    parser.add_argument("directory", type=str, help="削除する監視ディレクトリ")


def handle(args: argparse.Namespace):
    if args.command != "remove":
        return
    remove_scope(args.directory)
