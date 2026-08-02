"""Split selected residues into separate chains."""

import argparse
import string

import MDAnalysis as mda


def main() -> None:
    """Run the split2chain command."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_structure",
        type=str,
        help="Input structure file. e.g. PDB, CIF, etc.",
    )
    parser.add_argument(
        "residue_ids",
        type=str,
        help=(
            "Residue selections for each chain to extract. e.g. "
            '"1-20 30-40, 60-70". chain A: 1-20, 30-40; chain B: 60-70.'
        ),
    )

    args = parser.parse_args()

    u = mda.Universe(args.input_structure)

    chains = [
        u.select_atoms(f"protein and resid {res_ids}")
        for res_ids in args.residue_ids.split(",")
    ]

    for chain, chain_id in zip(chains, string.ascii_uppercase, strict=False):
        for atom in chain.atoms:
            atom.chainID = chain_id

    new_u = mda.Merge(*chains)
    if new_u.atoms is None:
        msg = "No atoms found in merged chains."
        raise ValueError(msg)
    new_u.atoms.write(f"{args.input_structure[:-4]}_merged.pdb")


if __name__ == "__main__":
    main()
