import argparse
from dit.model.extension import Extension


def _add(args):
    extension = Extension()

    for ext in args.extension:
        extension.add(ext)


def _remove(args):
    extension = Extension()

    for ext in args.extension:
        if not ext in extension.extensions:
            print(f"Error: 拡張子 {ext} は登録されていません")
            print("登録されているパターンは下の通りです")
            for registered_ext in extension.extensions:
                print(registered_ext)
            continue

        extension.remove(ext)


def _list(args):
    extension = Extension()

    for ext in extension.extensions:
        print(ext)


def register_subparser(root_parser):
    parser = root_parser.add_parser("ext")
    subparser = parser.add_subparsers(required=True)

    add_parser = subparser.add_parser("add", help="新たに拡張子を追加する")
    add_parser.add_argument("extension", type=str, help="新たに加える拡張子", nargs="+")
    add_parser.set_defaults(func=_add)

    remove_parser = subparser.add_parser("remove", help="拡張子を削除する")
    remove_parser.add_argument("extension", type=str, help="削除する拡張子", nargs="+")
    remove_parser.set_defaults(func=_remove)

    list_parser = subparser.add_parser("list", help="登録されている拡張子を表示する")
    list_parser.set_defaults(func=_list)
