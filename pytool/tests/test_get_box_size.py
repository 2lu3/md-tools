from pytool.structure import get_box_size
import MDAnalysis as mda


def test_get_box_size():
    u = mda.Universe('./tests/data/6rks.pdb')
    boxsize = get_box_size(u)

    print(boxsize)

