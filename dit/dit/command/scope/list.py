import argparse
from dit.model.scope import Scope
from natsort import natsorted


def list_scope():
    scope = Scope()

    directories = scope.directories

    for directory in natsorted(directories):
        print(directory)


def register_subparser(subparser):
    parser = subparser.add_parser(
        "list", help="Scopeに登録されているディレクトリを表示する"
    )

def handle(args: argparse.Namespace):
    list_scope()
