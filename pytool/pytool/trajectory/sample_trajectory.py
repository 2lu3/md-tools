import MDAnalysis as mda
from MDAnalysis.analysis import align
import click
from loguru import logger

import warnings
warnings.filterwarnings("ignore")


def _save_frame(average_structure: align.AverageStructure, start_frame: int, end_frame: int, output_name):
    # Use the AverageStructure class to compute the average
    av = average_structure.run(
        start=start_frame, stop=end_frame, step=1
    )

    assert av.universe.atoms is not None
    # Save the average structure
    with mda.Writer(output_name, av.universe.atoms.n_atoms) as W:
        W.write(av.results.universe.atoms)

    logger.info(f"Saved frame {start_frame}-{end_frame} -> {output_name}")


def sample_trajectory(
        dcd_path: str, pdb_path: str, output_name: str ="structure", output_num: int=3, window_size: int=10
):
    u = mda.Universe(pdb_path, dcd_path)
    average_structure = align.AverageStructure(u, reference=u, ref_frame=0)
    trajectory_num = len(u.trajectory)

    if output_num == 1:
        _save_frame(
            average_structure,
            trajectory_num - window_size,
            trajectory_num - 1,
            output_name + ".pdb",
        )
        return
    elif output_num == 2:
        _save_frame(
            average_structure,
            0,
            window_size - 1,
            output_name + "_0.pdb",
        )
        _save_frame(
            average_structure,
            trajectory_num - window_size,
            trajectory_num - 1,
            output_name + "_1.pdb",
        )
        return

    # 等間隔にフレームを抽出
    _save_frame(average_structure, 0, window_size - 1, output_name + "_0.pdb")
    for i in range(output_num - 2):
        total_extracted_frames = window_size * (output_num - 2)
        total_unextracted_frames = trajectory_num - total_extracted_frames
        segment_length = total_unextracted_frames // (output_num - 1)
        index = segment_length * (i + 1) + window_size * i
        _save_frame(
            average_structure,
            index,
            index + window_size - 1,
            output_name + "_" + str(i + 1) + ".pdb",
        )
    _save_frame(
        average_structure,
        trajectory_num - window_size,
        trajectory_num - 1,
        output_name + "_" + str(output_num - 1) + ".pdb",
    )


@click.command()
@click.argument("pdb_path")
@click.argument("dcd_path")
@click.option("--output", "-o", default="structure", help="file name of output file")
@click.option("--output_num", "-n", default=3, help="number of output file")
@click.option("--average_num", "-a", default=10, help="number of average frame")
def command(
    pdb_path: str, dcd_path: str, output: str, output_num: int, average_num: int
):
    sample_trajectory(dcd_path, pdb_path, output, output_num, average_num)
