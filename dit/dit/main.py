import argparse
from .command import add, clean, scope


def main():
    parser = argparse.ArgumentParser("DitはGitとDVCを併用するためのツールです。")
    subparser = parser.add_subparsers()

    #add.register_subparser(subparser)
    #clean.register_subparser(subparser)
    scope.register_subparser(subparser)

    args = parser.parse_args()
    args.func(args)
    #add.handle(args)
    #clean.handle(args)


if __name__ == "__main__":
    main()
