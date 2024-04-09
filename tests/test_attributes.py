import molecularnodes as mn
import pytest
import itertools


from .utils import sample_attribute_to_string
from .constants import (
    codes,
    attributes,
    data_dir
)

mn.unregister()
mn.register()

formats = ['pdb', 'cif', 'bcif']


@pytest.mark.parametrize("code, format", itertools.product(codes, formats))
def test_attribute(snapshot, code, format):
    mol = mn.io.fetch(code, cache_dir=data_dir, style=None, format=format)
    for attribute in attributes:
        snapshot.assert_match(
            sample_attribute_to_string(mol, attribute),
            f"att_{attribute}_values.txt"
        )
