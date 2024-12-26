import argparse
from .command import add, clean, scope

def main():
    parser = argparse.ArgumentParser("DitはGitとDVCを併用するためのツールです。")
    subparser = parser.add_subparsers(dest="command", required=True)

    add.register_subparser(subparser)
    clean.register_subparser(subparser)
    scope.register_subparser(subparser)

    args = parser.parse_args()
    if args.command == "add":
        add.handle(args)
    elif args.command == "clean":
        clean.handle(args)
    elif args.command == "scope":
        scope.handle(args)



if __name__ == "__main__":
    main()


