import MDAnalysis as mda
from pytool.structure import get_box_size


def test_get_box_size():
    u = mda.Universe("./tests/data/6rks.pdb")
    get_box_size(u)


