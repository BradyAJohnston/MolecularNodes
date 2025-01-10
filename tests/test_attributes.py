import itertools

import numpy as np
import pytest

import molecularnodes as mn

from .constants import attributes, codes, data_dir


formats = ["pdb", "cif", "bcif"]


@pytest.mark.parametrize("code, format", itertools.product(codes, formats))
def test_attribute(snapshot_custom, code, format):
    mol = mn.entities.fetch(code, cache_dir=data_dir, style=None, format=format)
    for attribute in attributes:
        try:
            assert snapshot_custom == mol.named_attribute(attribute)
        except AttributeError as e:
            assert snapshot_custom == e


def test_store_named_attribute(snapshot_custom):
    mol = mn.entities.fetch("8H1B", cache_dir=data_dir, style=None, format="bcif")
    before = mol.named_attribute("position")
    mol.store_named_attribute(mol.named_attribute("position") + 10, "position")
    after = mol.named_attribute("position")

    assert not np.allclose(before, after)
