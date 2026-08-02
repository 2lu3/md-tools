"""Reduce DCD trajectory frame counts."""

import argparse

from MDAnalysis.lib.formats.libdcd import DCDFile


def reduce_dcd(input_path: str, stride: int, out_path: str) -> None:
    """Write every stride-th frame from a DCD file.

    Args:
        input_path: Input DCD file path.
        stride: Frame stride.
        out_path: Output DCD file path.
    """
    with DCDFile(input_path) as dcd, DCDFile(
        out_path,
        "w",
    ) as out:
        out.write_header(
            remarks=dcd.header["remarks"],
            natoms=dcd.header["natoms"],
            istart=dcd.header["istart"] // stride,
            nsavc=dcd.header["nsavc"] // stride,
            delta=dcd.header["delta"] * stride,
            is_periodic=dcd.header["is_periodic"],
        )

        for i, ts in enumerate(dcd):
            if i % stride == 0:
                out.write(ts[0], ts[1])


def reduce_dcd_to_command() -> None:
    """Run the reduce-dcd command."""
    parser = argparse.ArgumentParser(
        description="Reduce the number of frames in a DCD file"
    )
    parser.add_argument("dcd", help="DCD file to reduce")
    parser.add_argument("stride", type=int, help="Stride to use")
    parser.add_argument("out", help="Output DCD file")

    args = parser.parse_args()

    reduce_dcd(args.dcd, args.stride, args.out)
