from argparse import ArgumentParser, Namespace
import subprocess
def commit(message: str):
    result = subprocess.run(["git", "commit", "-m", message], capture_output=True, text=True)
    print(result.stdout)

def register_subparser(subparser):
    parser = subparser.add_parser("commit", help="git commitを行う。")

    parser.add_argument("message", help="メッセージ", type=str)


def handle(args: Namespace):
    commit(args.message)

