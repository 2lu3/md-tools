import argparse
from .command import add

def main():
    parser = argparse.ArgumentParser("DitはGitとDVCを併用するためのツールです。")
    subparser = parser.add_subparsers(dest="command", required=True)

    add.register_subparser(subparser)

    args = parser.parse_args()
    if args.command == "add":
        add.handle(args)


if __name__ == "__main__":
    main()


