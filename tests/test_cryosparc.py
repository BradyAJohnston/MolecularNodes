import numpy as np
import pytest
from molecularnodes.entities.ensemble.cryosparc import Dataset
from .constants import data_dir


@pytest.fixture
def dset():
    d = np.load(data_dir / "cryosparc" / "J123_particles_exported.cs")
    return Dataset(dset=d)


class TestDataset:
    def test_cryosparc_dset_len(self, dset):
        assert len(dset) == 100

    def test_cryosparc_dset_fields(self, dset):
        d = np.load(data_dir / "cryosparc" / "J123_particles_exported.cs")
        for field_name in d.dtype.names:
            assert field_name in dset
            assert np.array_equal(d[field_name], dset[field_name])
            assert np.array_equal(dset[field_name], dset.get(field_name))

        assert dset.get("BAD_FIELD") is None
        assert np.array_equal(dset.get("BAD_FIELD", np.arange(3)), np.arange(3))

    def test_cryosparc_dset_bob_entry(self, dset):
        field_name = "location/center_x_frac"
        field_type = "FLOAT"
        bob = dset.bob_entry(field_name, field_type)
        assert bob["name"] == field_name
        assert np.array_equal(bob["data"], dset[field_name])
        assert bob["atype"] == field_type
