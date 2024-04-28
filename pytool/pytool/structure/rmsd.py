import click
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import matplotlib.pyplot as plt

import warnings

warnings.filterwarnings("ignore", "DCDReader currently makes independent timesteps")


def rmsd(u: mda.Universe, ref: mda.Universe, select: str = "protein") -> float:
    u_selected = u.select_atoms(select)
    ref_selected = ref.select_atoms(select)

    return rms.rmsd(
        u_selected.positions,
        ref_selected.positions,
        center=True,
        superposition=True,
    )


def rmsd_trajectory(
    u: mda.Universe, ref: mda.Universe, select="protein"
) -> tuple[list[float], list[float]]:
    R = rms.RMSD(
        u,
        ref,
        select=select,
        verbose=True,
    )

    R.run()

    result = R.results.rmsd.T
    time = result[1]
    rmsd = result[2]
    return time, rmsd


@click.command()
@click.argument("toplogy1", type=click.Path(exists=True))
@click.argument("toplogy2", type=click.Path(exists=True))
@click.option("--select", default="protein", help="RMSDを計算する構造を指定する")
def rmsd_to_command(toplogy1: str, toplogy2: str, select: str):
    u = mda.Universe(toplogy1)
    ref = mda.Universe(toplogy2)

    print(rmsd(u, ref, select=select))


@click.command()
@click.argument("toplogy1", type=click.Path(exists=True))
@click.argument("toplogy2", type=click.Path(exists=True))
@click.option("--select", default="protein", help="RMSDを計算する構造を指定する")
@click.option("--picture", default="rmsd.png", help="出力画像ファイル名")
def rmsd_trajectory_to_command(toplogy1: str, toplogy2: str, select: str, picture: str):
    u = mda.Universe(toplogy1)
    ref = mda.Universe(toplogy2)

    time, rmsd = rmsd_trajectory(u, ref, select=select)

    fig: plt.Figure = plt.figure()  # type: ignore
    ax = fig.add_subplot(111)

    ax.plot(time, rmsd, label="all")
    ax.legend(loc="best")
    ax.set_xlabel("time (ps)")
    ax.set_ylabel(r"RMSD ($\AA$)")
    fig.savefig(picture)
