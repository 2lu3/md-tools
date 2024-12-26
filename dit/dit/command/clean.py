import argparse

# dic gc -w --not-in-remote

def register_subparser(subparser):
    parser = subparser.add_parser("clean", help="ストレージを可能な限り空ける")

def handle(args: argparse.Namespace):
    if args.command != "clean":
        return
    pass
