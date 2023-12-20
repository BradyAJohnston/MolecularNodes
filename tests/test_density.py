import molecularnodes as mn
import pytest
import numpy as np
from .constants import test_data_directory
try:
    import pyopenvdb
except ImportError:
    pytest.skip("pyopenvdb not installed", allow_module_level=True)

mn.unregister()
mn.register()

def test_density_load():
    file = test_data_directory / "emd_24805.map.gz"
    obj = mn.io.density.load(file,style="density_surface")
    assert obj.mn.molecule_type == "density"