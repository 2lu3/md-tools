import click
from .command.scope import scope
from .command.extension import ext
from .command.add import add
from .command.commit import commit
from .command.push import push


@click.group()
def cli():
    pass


def main():
    cli.add_command(scope)
    cli.add_command(ext)
    cli.add_command(add)
    cli.add_command(commit)
    cli.add_command(push)

    cli()


if __name__ == "__main__":
    main()
