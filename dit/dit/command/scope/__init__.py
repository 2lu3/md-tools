import argparse

def register_subparser(subparser):
    parser = subparser.add_parser("scope", help="Scopeに関する操作")

    subsubparser = parser.add_subparsers(dest="command", help="Scopeに関する操作")

    from .add import register_subparser as add_register_subparser
    add_register_subparser(subsubparser)

    from .list import register_subparser as list_register_subparser
    list_register_subparser(subsubparser)

    from .remove import register_subparser as remove_register_subparser
    remove_register_subparser(subsubparser)



def handle(args: argparse.Namespace):
    if args.command != "scope":
        return

    from .add import handle as add_handle
    from .list import handle as list_handle
    from .remove import handle as remove_handle

    add_handle(args)
    list_handle(args)
    remove_handle(args)

