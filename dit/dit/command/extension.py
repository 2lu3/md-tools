import click
from dit.model.extension import Extension


@click.group()
def ext():
    pass


@ext.command()
@click.argument("extension")
def add(extension):
    Extension().add(extension)

@ext.command()
@click.argument("extension")
def remove(extension):
    ext_manager = Extension()

    if not extension in ext_manager.extensions:
        print(f"Error: 拡張子 {ext} は登録されていません")
        print("登録されているパターンは下の通りです")
        for registered_ext in ext_manager.extensions:
            print(registered_ext)
        return

    ext_manager.remove(extension)


@ext.command()
def list():
    extension = Extension()

    for ext in extension.extensions:
        print(ext)
