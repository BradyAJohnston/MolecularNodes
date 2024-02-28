import molecularnodes as mn
import pytest
import itertools
import tempfile


from .utils import sample_attribute_to_string
from .constants import (
    codes,
    attributes
)

mn.unregister()
mn.register()

formats = ['pdb', 'mmtf', 'cif']

with tempfile.TemporaryDirectory() as temp_dir:
    @pytest.mark.parametrize("code, format", itertools.product(codes, formats))
    def test_attribute(snapshot, code, format):
        mol = mn.io.fetch(code, cache_dir=temp_dir, style=None, format=format)
        for attribute in attributes:
            snapshot.assert_match(
                sample_attribute_to_string(mol, attribute),
                f"att_{attribute}_values.txt"
            )
