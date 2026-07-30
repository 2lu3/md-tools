"""Set chain identifiers in PDB files."""

from pathlib import Path

import click


def set_chain_id(pdb_path: str, output_path: str, chain_id: str) -> None:
    """Set the chain identifier for ATOM records in a PDB file.

    Args:
        pdb_path: Input PDB file path.
        output_path: Output PDB file path.
        chain_id: Chain identifier to write.
    """
    lines = []
    with Path(pdb_path).open() as f:
        for line in f:
            if line.startswith("ATOM"):
                lines.append(line[:21] + chain_id + line[22:])
            else:
                lines.append(line)

    with Path(output_path).open("w") as f:
        f.writelines(lines)


@click.command()
@click.option("--pdb_path", "-i", help="Input pdb file path")
@click.option("--output_path", "-o", help="Output pdb file path")
@click.option("--chain_id", "-c", help="Chain ID")
def set_chain_id_to_command(pdb_path: str, output_path: str, chain_id: str) -> None:
    """Run the set-chain-id command."""
    set_chain_id(pdb_path, output_path, chain_id)
