import argparse
from MDAnalysis.lib.formats.libdcd import DCDFile


def reduce_dcd(input_path: str, stride: int, out_path: str):
    with DCDFile(input_path) as dcd:
        print(dcd.header)
        with DCDFile(
            out_path,
            "w",
        ) as out:
            print(dcd.header)
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


def main():
    parser = argparse.ArgumentParser(
        description="Reduce the number of frames in a DCD file"
    )
    parser.add_argument("dcd", help="DCD file to reduce")
    parser.add_argument("stride", type=int, help="Stride to use")
    parser.add_argument("out", help="Output DCD file")

    args = parser.parse_args()

    reduce_dcd(args.dcd, args.stride, args.out)
