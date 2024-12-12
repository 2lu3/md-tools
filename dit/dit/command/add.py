import argparse
import os

def add(path: list[str], is_all: bool):
    pass

def to_command():
    parser = argparse.ArgumentParser()

    parser.add_argument("path", help="Path to the file or directory to add.", type=str, nargs='*')
    parser.add_argument("-A", "--all", help="Add all files in the scope.", action="store_true")

    args = parser.parse_args()

    add(args.path, args.all)

