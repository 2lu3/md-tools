from typing import Optional, Union
import MDAnalysis as mda
from MDAnalysis.analysis.align import AlignTraj, alignto
import click

import warnings

warnings.filterwarnings("ignore")


def align(
    mobile: mda.Universe,
    ref: mda.Universe,
    output_path: str,
    select: str = "protein and name CA",
):
    alignto(mobile, ref, select=select, weights="mass")

    output = mda.Writer(output_path, mobile.atoms.n_atoms)  # type: ignore
    output.write(mobile.atoms)


@click.command()
@click.argument("target", type=click.Path(exists=True))
@click.argument("ref", type=click.Path(exists=True))
@click.option("--output_path", "-o", help="Path to output file", required=True)
@click.option("--select", "-s", help="Selection string", default="protein and name CA")
def align_to_command(
    target: str, ref: str, output_path: str, select: str = "protein and name CA"
):
    align(mda.Universe(target), mda.Universe(ref), output_path, select)


def align_trajectory(
    mobile: mda.Universe,
    ref: mda.Universe,
    output_path: str,
    select: str = "protein and name CA",
    verbose: bool = True,
):
    AlignTraj(mobile, ref, select=select, filename=output_path, verbose=verbose).run()


@click.command()
@click.option("--pdb_path", "-p", help="Path to PDB file", required=True)
@click.option("--dcd_path", "-d", help="Path to DCD file", required=True)
@click.option("--output_path", "-o", help="Path to output file", required=True)
@click.option("--select", "-s", help="Selection string", default="protein and name CA")
def align_trajectory_to_command(
    pdb_path: str, dcd_path: str, output_path: str, select: str = "protein and name CA"
):
    align_trajectory(
        mda.Universe(pdb_path, dcd_path), mda.Universe(pdb_path), output_path, select
    )
