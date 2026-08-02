import MDAnalysis as mda
from pytool.sequence import rename2charmm


def test_rename2charmm():
    u = mda.Universe("./tests/data/6rks.pdb")
    rename2charmm(u)
