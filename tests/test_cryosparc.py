import numpy as np
import pytest
from molecularnodes.entities.ensemble.cryosparc import (
    MNDataset,
    i32_vec_to_uid,
    uid_as_i32_vec,
)
from .constants import data_dir

CS_FILE_NAME = "J123_particles_exported.cs"


@pytest.fixture
def darray():
    return np.load(data_dir / "cryosparc" / CS_FILE_NAME)


@pytest.fixture
def dset():
    d = np.load(data_dir / "cryosparc" / CS_FILE_NAME)
    return MNDataset(dset=d)


def drop_fields(a: np.typing.NDArray, *fields) -> np.typing.NDArray:
    if a.dtype.names is None:
        raise AttributeError("Array is not structured")
    for f in fields:
        if f not in a.dtype.names:
            raise ValueError(f"{f} not in array names")
    return a[[f for f in a.dtype.names if f not in fields]]


class TestDataset:
    def test_dset_len(self, dset, darray):
        assert len(dset) == len(darray)

    def test_dset_fields(self, darray, dset):
        for field_name in darray.dtype.names:
            assert field_name in dset
            assert np.array_equal(darray[field_name], dset[field_name])
            assert np.array_equal(dset[field_name], dset.get(field_name))

        assert dset.get("BAD_FIELD") is None
        assert np.array_equal(dset.get("BAD_FIELD", np.arange(3)), np.arange(3))

    def test_dset_bob_entry(self, dset):
        field_name = "location/center_x_frac"
        field_type = "FLOAT"
        bob = dset.bob_entry(field_name, field_type)
        assert bob["name"] == field_name
        assert np.array_equal(bob["data"], dset[field_name])
        assert bob["atype"] == field_type

    def test_uid_as_i32(self, dset):
        split_uids = uid_as_i32_vec(dset["uid"])
        assert split_uids.shape == (len(dset), 2)
        combined_uids = i32_vec_to_uid(split_uids=split_uids)
        assert np.array_equal(combined_uids, dset["uid"])

    def test_shift2d(self, darray, dset):
        assert np.allclose(
            dset.shift2d,
            darray["alignments2D/shift"]
            * darray["alignments2D/psize_A"][:, np.newaxis],
        )
        dset.dset = drop_fields(dset.dset, "alignments2D/psize_A")
        assert dset.shift2d is None

    def test_shift3d(self, darray, dset):
        assert np.allclose(
            dset.shift3d,
            darray["alignments3D/shift"]
            * darray["alignments3D/psize_A"][:, np.newaxis],
        )
        dset.dset = drop_fields(dset.dset, "alignments3D/psize_A")
        assert dset.shift3d is None

    def test_positions(self, dset, darray):
        positions = dset.position
        mic_size = (
            darray["location/micrograph_shape"]
            * darray["location/micrograph_psize_A"][:, np.newaxis]
        )
        darray_positions = mic_size * np.column_stack(
            (darray["location/center_x_frac"], darray["location/center_y_frac"])
        )
        darray_defocus = np.mean([darray["ctf/df1_A"], darray["ctf/df2_A"]], axis=0)
        darray_defocus = np.median(darray_defocus) - darray_defocus
        assert np.allclose(positions[:, :2], darray_positions)
        assert np.allclose(positions[:, 2], darray_defocus)

        dset.dset = drop_fields(dset.dset, "location/center_x_frac")
        assert np.array_equal(
            dset.position[:, 0], np.zeros(len(dset), dtype=np.float32)
        )
        assert np.array_equal(dset.position[:, 1], positions[:, 1])

        dset.dset = drop_fields(dset.dset, "location/center_y_frac")
        assert np.allclose(dset.position[:, :2], 0)

        dset.dset = darray.copy()
        dset.dset = drop_fields(dset.dset, "location/center_y_frac")
        assert np.array_equal(dset.position[:, 0], positions[:, 0])
        assert np.array_equal(
            dset.position[:, 1], np.zeros(len(dset), dtype=np.float32)
        )

        dset.dset = drop_fields(dset.dset, "ctf/df1_A", "ctf/df2_A")
        assert np.allclose(dset.position[:, 2], 0)
