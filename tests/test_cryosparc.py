import re
import numpy as np
import pytest
from cryosparc.dataset import Dataset
from molecularnodes.entities.ensemble.cryosparc import (
    UID_REGEX,
    MNDataset,
    i32_vec_to_uid,
    uid_as_i32_vec,
)
from .constants import data_dir

CS_FILE_NAME = "J123_particles_exported.cs"


@pytest.fixture
def tools_dset():
    return Dataset.load(data_dir / "cryosparc" / CS_FILE_NAME)


@pytest.fixture
def dset():
    d = Dataset.load(data_dir / "cryosparc" / CS_FILE_NAME)
    return MNDataset(dset=d)


class TestDataset:
    def test_dset_len(self, dset, tools_dset):
        assert len(dset) == len(tools_dset)

    def test_dset_fields(self, tools_dset, dset):
        for field_name in tools_dset.fields():
            assert field_name in dset
            assert np.array_equal(tools_dset[field_name], dset[field_name])
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

    def test_shift2d(self, tools_dset, dset):
        assert np.allclose(
            dset.shift2d,
            tools_dset["alignments2D/shift"]
            * tools_dset["alignments2D/psize_A"][:, np.newaxis],
        )
        dset.dset.drop_fields(["alignments2D/psize_A"])
        assert dset.shift2d is None

    def test_shift3d(self, tools_dset, dset):
        assert np.allclose(
            dset.shift3d,
            tools_dset["alignments3D/shift"]
            * tools_dset["alignments3D/psize_A"][:, np.newaxis],
        )
        dset.dset.drop_fields(["alignments3D/psize_A"])
        assert dset.shift3d is None

    def test_positions(self, dset, tools_dset):
        positions = dset.position
        mic_size = (
            tools_dset["location/micrograph_shape"]
            * tools_dset["location/micrograph_psize_A"][:, np.newaxis]
        )
        tools_dset_positions = mic_size * np.column_stack(
            (tools_dset["location/center_x_frac"], tools_dset["location/center_y_frac"])
        )
        tools_dset_defocus = np.mean(
            [tools_dset["ctf/df1_A"], tools_dset["ctf/df2_A"]], axis=0
        )
        tools_dset_defocus = np.median(tools_dset_defocus) - tools_dset_defocus
        assert np.allclose(positions[:, :2], tools_dset_positions)
        assert np.allclose(positions[:, 2], tools_dset_defocus)

        dset.dset.drop_fields(["location/center_x_frac"])
        assert np.array_equal(
            dset.position[:, 0], np.zeros(len(dset), dtype=np.float32)
        )
        assert np.array_equal(dset.position[:, 1], positions[:, 1])

        dset.dset.drop_fields(["location/center_y_frac"])
        assert np.allclose(dset.position[:, :2], 0)

        dset.dset = tools_dset.copy()
        dset.dset.drop_fields(["location/center_y_frac"])
        assert np.array_equal(dset.position[:, 0], positions[:, 0])
        assert np.array_equal(
            dset.position[:, 1], np.zeros(len(dset), dtype=np.float32)
        )

        dset.dset.drop_fields(["ctf/df1_A", "ctf/df2_A"])
        assert np.allclose(dset.position[:, 2], 0)

    def test_blob_uids(self, dset):
        cs_blob_uids = np.array(
            # ignore type b/c if the test dataset is missing UIDs we want to know about it
            [re.search(UID_REGEX, s).group(1) for s in dset.dset["blob/path"]],  # type: ignore
            dtype=np.uint64,
        )
        assert np.array_equal(i32_vec_to_uid(dset.blob_path_uids), cs_blob_uids)
        dset.dset["blob/path"][0] = "test missing uid"
        assert dset.blob_path_uids is None
        dset.dset.drop_fields(["blob/path"])
        assert dset.blob_path_uids is None
