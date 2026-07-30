"""Pick residue mappings from aligned FASTA sequences."""

from argparse import ArgumentParser
from dataclasses import dataclass

from Bio import SeqIO


@dataclass
class Residue:
    """A residue in an aligned sequence."""

    id: int | None
    name: str


class Sequence:
    """An aligned sequence with residue identifiers."""

    def __init__(self, sequence: str, offset: int) -> None:
        """Initialize a sequence with a residue ID offset.

        Args:
            sequence: Aligned sequence text.
            offset: Residue ID offset for non-gap residues.
        """
        self.sequence = sequence
        self.offset = offset
        self.amino_acids = self._create_amino_acids()

    def _create_amino_acids(self) -> list[Residue]:
        amino_acids = []
        index = 0
        for residue_name in self.sequence:
            if residue_name == ".":
                amino_acids.append(Residue(None, residue_name))
            else:
                amino_acids.append(Residue(index + self.offset, residue_name))
                index += 1
        return amino_acids

    def __str__(self) -> str:
        """Return a compact sequence summary."""
        dot_count = 0
        for amino_acid in self.amino_acids:
            if amino_acid.name == ".":
                dot_count += 1
        return f"Sequence<dot={dot_count} amino={len(self.amino_acids) - dot_count}>"


def load_fasta(path: str, ref_index: int, target_index: int) -> tuple[str, str]:
    """Load reference and target sequences from a FASTA file.

    Args:
        path: FASTA file path.
        ref_index: Reference sequence index.
        target_index: Target sequence index.

    Returns:
        Reference and target sequences.
    """
    records = list(SeqIO.parse(path, "fasta"))
    return str(records[ref_index].seq), str(records[target_index].seq)


def main() -> None:
    """Run the pickup-aligned-seq command."""
    parser = ArgumentParser()
    parser.add_argument("input", type=str, help="FASTA file")
    parser.add_argument("--ref-index", type=int, help="Reference sequence index.")
    parser.add_argument("--target-index", type=int, help="Target sequence index.")
    parser.add_argument(
        "--target-offset",
        type=int,
        help="Residue ID of the first target amino acid.",
    )
    parser.add_argument(
        "--ref-offset",
        type=int,
        help="Residue ID of the first reference amino acid.",
    )
    parser.add_argument(
        "--ref-residue-ids",
        type=str,
        help=(
            "Comma-separated reference residue IDs corresponding to target "
            "sequence residues. e.g. 1,2,3."
        ),
    )

    args = parser.parse_args()

    ref_seq_str, target_seq_str = load_fasta(
        args.input, args.ref_index, args.target_index
    )

    ref_sequence = Sequence(ref_seq_str, args.ref_offset)
    target_sequence = Sequence(target_seq_str, args.target_offset)

    ref_searching_residue_ids = [int(x) for x in args.ref_residue_ids.split(",")]

    result: dict[str, str] = {}
    for i in range(len(ref_sequence.amino_acids)):
        if ref_sequence.amino_acids[i].id in ref_searching_residue_ids:
            result[
                str(ref_sequence.amino_acids[i].id)
            ] = (
                f"{ref_sequence.amino_acids[i].id} "
                f"{ref_sequence.amino_acids[i].name} "
                f"{target_sequence.amino_acids[i].id} "
                f"{target_sequence.amino_acids[i].name}"
            )

    for _ref_residue_id in ref_searching_residue_ids:
        pass
