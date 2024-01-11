import molecularnodes as mn
import pytest
from .utils import sample_attribute_to_string
from .constants import (
    codes,
    attributes
)


@pytest.mark.parametrize("code", codes)
def test_attribute(snapshot, code, tmpdir):
    mol = mn.io.pdb.load(code, cache_dir=tmpdir)
    for attribute in attributes:
        snapshot.assert_match(
            sample_attribute_to_string(mol, attribute),
            f"att_{attribute}_values.txt"
        )
