import molecularnodes as mn
import pytest
import numpy as np
from .utils import sample_attribute_to_string
from .constants import test_data_directory

def test_starfile_attributes(snapshot):
    for type in ["cistem", "relion"]:
        file = test_data_directory / f"{type}.star"
        obj = mn.io.star.load(file)
        for attribute in obj.data.attributes.keys():
            snapshot.assert_match(
                sample_attribute_to_string(obj, attribute, n = 200), 
                f"{type}_{attribute}_values"
            )