from typing import Optional
import matplotlib.pyplot as plt
import click

def plot_2d(filename: str, x_index:int, y_indexes: list[int], title: str, xlabel: str, ylabel: str, save_path: Optional[str]= None):
    with open(filename, 'r') as f:
        lines = f.readlines()
        x = []
        ys = [[] for _ in range(len(y_indexes))]
        for line in lines:
            line = line.strip()
            if line.startswith('#'):
                continue
            values = line.split()
            x.append(float(values[x_index]))
            for i, y_index in enumerate(y_indexes):
                ys[i].append(float(values[y_index]))
        for i, y in enumerate(ys):
            plt.plot(x, y, label=f'column {y_indexes[i]}')
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.legend()

        if save_path is not None:
            plt.savefig(save_path)
        else:
            plt.show()

@click.command()
@click.argument('filename')
@click.argument('x_index', type=int)
@click.argument('y_indexes', nargs=-1, type=int)
@click.option('--title', '-t', default='Plot', help='Title of the plot')
@click.option('--xlabel', '-x', default='x', help='Label of x axis')
@click.option('--ylabel', '-y', default='y', help='Label of y axis')
@click.option('--output_path', '-o', default=None, help='Path to save the plot')
def plot_2d_to_command(filename: str, x_index:int, y_indexes: list[int], title: str, xlabel: str, ylabel: str, output_path: Optional[str]):
    plot_2d(filename, x_index, y_indexes, title, xlabel, ylabel, output_path)
