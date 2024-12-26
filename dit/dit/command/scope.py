import argparse
from dit.model.scope import Scope
import os
from natsort import natsorted


def add(args):
    scope = Scope()
    assert os.path.isdir(args.directory), f"{args.directory} is not a directory"
    scope.add_directory(args.directory)


def list():
    scope = Scope()

    directories = scope.directories

    for directory in natsorted(directories):
        print(directory)

def remove(args):
    scope = Scope()
    scope.remove_directory(args.directory)

def register_subparser(subparser):
    parser = subparser.add_subparsers("scope", help="scope subcommands")

    subsubparser = parser.add_subparsers(dest="command", help="scope subcommands")

    add_parser = subsubparser.add_parser("add", help="add directory to scope")
    add_parser.add_argument("directory", help="directory to add")
    add_parser.set_defaults(handler=add)

    list_parser = subsubparser.add_parser("list", help="list directories in scope")
    list_parser.set_defaults(handler=list)

    remove_parser = subsubparser.add_parser("remove", help="remove directory from scope")
    remove_parser.add_argument("directory", help="directory to remove")
    remove_parser.set_defaults(handler=remove)

