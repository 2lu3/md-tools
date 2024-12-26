import argparse
from dit.model.scope import Scope
import os
from natsort import natsorted
from loguru import logger


def add_scope(args):
    scope = Scope()
    assert os.path.isdir(args.directory), f"{args.directory} is not a directory"
    scope.add_directory(args.directory)
    logger.info(f"add {args.directory} to scope")


def list_scope(args):
    scope = Scope()

    directories = scope.directories

    for directory in natsorted(directories):
        print(directory)

def remove_scope(args):
    scope = Scope()
    scope.remove_directory(args.directory)

def register_subparser(subparser):
    parser = subparser.add_parser("scope")
    subparser = parser.add_subparsers()


    add_parser = subparser.add_parser("add", help="add directory to scope")
    add_parser.add_argument("directory", help="directory to add")
    add_parser.set_defaults(func=add_scope)

    list_parser = subparser.add_parser("list", help="list directories in scope")
    list_parser.set_defaults(func=list_scope)

    remove_parser = subparser.add_parser("remove", help="remove directory from scope")
    remove_parser.add_argument("directory", help="directory to remove")
    remove_parser.set_defaults(func=remove_scope)


