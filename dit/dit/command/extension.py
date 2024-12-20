import argparse
from dit.model.extension import Extension

def add_extension_to_command():
    parser = argparse.ArgumentParser()

    parser.add_argument("extension", type=str, help="新たに加える拡張子", nargs="+")

    args = parser.parse_args()

    extension = Extension()

    for extension in args.extension:
        extension.add(extension)


def remove_extension_to_command():
    parser = argparse.ArgumentParser()

    parser.add_argument("extension", type=str, help="削除する拡張子", nargs="+")

    args = parser.parse_args()

    extension = Extension()

    for extension in args.extension:
        if not extension in extension.extensions:
            print(f"Error: 拡張子 {extension} は登録されていません")
            print("登録されているパターンは下の通りです")
            for registered_ext in extension.extensions:
                print(registered_ext)
            continue

        extension.add_extension(extension)

def list_extension_to_command():
    extension = Extension()

    for ext in extension.extensions:
        print(ext)

def from_parser(parser: argparse.ArgumentParser):
    subparser = parser.add_subparsers()

    add_parser = subparser.add_parser("add", help="新たに拡張子を追加する")

    remove_parser = subparser.add_parser("remove", help="拡張子を削除する")

    list_parser = subparser.add_parser("list", help="登録されている拡張子を表示する")

