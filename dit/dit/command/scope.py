import argparse
from dit.model.scope import Scope
import os


def add(args):
    scope = Scope()
    assert os.path.isdir(args.directory), f"{args.directory} is not a directory"
    scope.add_directory(args.directory)


