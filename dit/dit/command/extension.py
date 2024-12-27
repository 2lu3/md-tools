import click
from dit.model.extension import Extension


@click.group()
def ext():
    pass


@ext.command()
@click.argument("extension")
def add(ext):
    Extension().add(ext)

@ext.command()
@click.argument("extension")
def remove(ext):
    extension = Extension()

    if not ext in extension.extensions:
        print(f"Error: 拡張子 {ext} は登録されていません")
        print("登録されているパターンは下の通りです")
        for registered_ext in extension.extensions:
            print(registered_ext)
        return

    extension.remove(ext)


@ext.command()
def list():
    extension = Extension()

    for ext in extension.extensions:
        print(ext)
