import click
from .command.scope import scope
from .command.extension import ext

@click.group()
def cli():
    pass


def main():
    cli.add_command(scope)
    cli.add_command(ext)

    cli()


if __name__ == "__main__":
    main()
