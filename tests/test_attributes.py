import molecularnodes as mn
import pytest
from .utils import sample_attribute_to_string
from .constants import (
    codes, 
    attributes
)

@pytest.mark.parametrize("code", codes)
def test_attribute(snapshot, code):
    mol =  mn.load.molecule_rcsb(code)
    for attribute in attributes:
        snapshot.assert_match(
            sample_attribute_to_string(mol, attribute), 
            f"att_{attribute}_values.txt"
        )