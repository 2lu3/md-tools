import click
from .command.scope import scope
from .command.extension import ext
from .command.add import add
from .command.commit import commit
from .command.push import push
from .command.pull import pull
from .command.clean import clean


@click.group()
def cli():
    pass


def main():
    cli.add_command(scope)
    cli.add_command(ext)
    cli.add_command(add)
    cli.add_command(commit)
    cli.add_command(push)
    cli.add_command(pull)
    cli.add_command(clean)

    cli()


if __name__ == "__main__":
    main()
