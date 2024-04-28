from pytool.sequence import rename2charmm
import MDAnalysis as mda

def test_rename2charmm():
    u = mda.Universe('./tests/data/6rks.pdb')
    new_u = rename2charmm(u)
