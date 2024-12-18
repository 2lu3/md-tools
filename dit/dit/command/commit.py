from argparse import ArgumentParser
def commit(message: str):
    pass

def to_command():
    parser = ArgumentParser()

    parser.add_argument("message", help="メッセージ", type=str)
    args = parser.parse_args()

    commit(args.message)

