from argparse import ArgumentParser
from typing import Optional
from Bio import SeqIO
from dataclasses import dataclass


@dataclass
class Residue:
    id: Optional[int]
    name: str


class Sequence:
    def __init__(self, sequence: str, offset: int):
        self.sequence = sequence
        self.offset = offset
        self.amino_acids = self._create_amino_acids()

    def _create_amino_acids(self):
        amino_acids = []
        index = 0
        for residue_name in self.sequence:
            if residue_name == ".":
                amino_acids.append(Residue(None, residue_name))
            else:
                amino_acids.append(Residue(index + self.offset, residue_name))
                index += 1
        return amino_acids

    def __str__(self):
        dot_count = 0
        for amino_acid in self.amino_acids:
            if amino_acid.name == ".":
                dot_count += 1
        return f"Sequence<dot={dot_count} amino={len(self.amino_acids) - dot_count}>"



def load_fasta(path: str, ref_index: int, target_index: int):
    records = list(SeqIO.parse(path, "fasta"))
    return records[ref_index].seq, records[target_index].seq


def main():
    parser = ArgumentParser()
    parser.add_argument("input", type=str, help="FASTA file")
    parser.add_argument("--ref-index", type=int, help="参照とする配列のindex.")
    parser.add_argument("--target-index", type=int,help="目標とする配列のindex.")
    parser.add_argument("--target-offset", type=int, help="最初のアミノ酸のresidue id")
    parser.add_argument("--ref-offset", type=int, help="最初のアミノ酸のresidue id")
    parser.add_argument(
        "--ref-residue-ids",
        type=str,
        help="求めたいtarget sequenceに対応するref sequenceのresidue idをコンマ区切りで指定. e.g. 1,2,3",
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
            ] = f"{ref_sequence.amino_acids[i].id} {ref_sequence.amino_acids[i].name} {target_sequence.amino_acids[i].id} {target_sequence.amino_acids[i].name}"

    for ref_residue_id in ref_searching_residue_ids:
        print(result[str(ref_residue_id)])
