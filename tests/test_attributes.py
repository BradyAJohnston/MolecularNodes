import molecularnodes as mn
import pytest
import itertools
import numpy as np


from .utils import sample_attribute
from .constants import codes, attributes, data_dir

mn._test_register()

formats = ["pdb", "cif", "bcif"]


@pytest.mark.parametrize("code, format", itertools.product(codes, formats))
def test_attribute(snapshot_custom, code, format):
    mol = mn.io.fetch(code, cache_dir=data_dir, style=None, format=format)
    for attribute in attributes:
        vals = sample_attribute(mol, attribute)
        assert snapshot_custom == vals


def test_set_attribute(snapshot_custom):
    mol = mn.io.fetch("8H1B", cache_dir=data_dir, style="surface", format="bcif")
    before = sample_attribute(mol, "position")
    mol.set_attribute(mol.get_attribute("position") + 10, "position")
    after = sample_attribute(mol, "position")

    assert not np.allclose(before, after)
