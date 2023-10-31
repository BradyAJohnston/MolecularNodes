import molecularnodes as mn
from .utils import sample_attribute_to_string
from .constants import (
    test_data_directory
)

def test_cellpack_data(snapshot):
    object, collection = mn.pack.open_file(
        test_data_directory / "synvesicle_2-no_bonds.bcif"
    )
    attributes = object.data.attributes.keys()
    for attribute in attributes:
        snapshot.assert_match(
            sample_attribute_to_string(object, attribute), 
            f"att_{attribute}_values.txt"
        )