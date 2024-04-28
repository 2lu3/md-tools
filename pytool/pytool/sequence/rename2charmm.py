import click
import MDAnalysis as mda

def rename2charmm(u: mda.Universe):
    """Rename atoms in a universe to CHARMM naming convention

    Args:
        u (mda.Universe): A universe


    Rename "PDB general atom name" to "CHARMM-specific atom name"
        HIS => HSD
        CD1 atom of ILE => CD
        C-terminal carboxyl oxygen O and OXT => OT1 and OT2
    """

    # Rename CD1 atom of ILE to CD
    u.select_atoms("resname ILE and name CD1").names = "CD"

    # Rename HIS to HSD
    u.select_atoms("resname HIS").residues.resnames = "HSD"

    # C terminal residue for every chain
    for seg in u.segments:
        last_resid = seg.residues[-1].resid
        u.select_atoms(f"resid {last_resid} and name O").names = "OT1"
        u.select_atoms(f"resid {last_resid} and name OXT").names = "OT2"

    return u


@click.command()
@click.argument("input_pdb", type=click.Path(exists=True))
@click.argument("output_pdb", type=click.Path(exists=False))
def main(input_pdb: str, output_pdb: str):
    """Rename atoms in a PDB file to CHARMM naming convention

    Args:
        input_pdb (str): Input PDB file name
        output_pdb (str): Output PDB file name
    """

    u = mda.Universe(input_pdb)

    renamed_u = rename2charmm(u)

    renamed_u.atoms.write(output_pdb)
