import json
import databpy
import numpy as np
import pytest
import starfile
from nodebpy.nodes.geometry import StoreNamedAttribute
from pandas import DataFrame
from scipy.spatial.transform import Rotation
import molecularnodes as mn
from tests.utils import GeometrySet
from .constants import data_dir


@pytest.mark.parametrize("type", ["cistem", "relion"])
def test_starfile_attributes(type, snapshot):
    """
    Test that our nodes correctly convert the starfile attribute columns to quaternions that matches the convention of scipy.spatial.transform.Rotation.from_euler.
    """
    file = data_dir / f"starfile/{type}.star"
    ensemble = mn.entities.ensemble.StarFile.load(file)

    star = starfile.read(file)

    if type == "relion":
        assert isinstance(star, dict)
        df: DataFrame = star["particles"].merge(star["optics"], on="rlnOpticsGroup")  # ty: ignore[unresolved-attribute]
        euler_angles = df[["rlnAngleRot", "rlnAngleTilt", "rlnAnglePsi"]].to_numpy()

    elif type == "cistem":
        assert isinstance(star, DataFrame)
        df = star
        euler_angles = df[
            ["cisTEMAnglePhi", "cisTEMAngleTheta", "cisTEMAnglePsi"]
        ].to_numpy()

    # Calculate Scipy rotation from the euler angles
    # Note: rot_from_euler = quats
    rotation_scipy = Rotation.from_euler(
        seq="ZYZ", angles=euler_angles, degrees=True
    ).inv()

    with ensemble.tree.reset() as (atoms, join):
        rot = {
            "cistem": mn.nodes.geometry.RotationCisTEM,
            "relion": mn.nodes.geometry.RotationRELION,
        }

        (
            atoms
            >> StoreNamedAttribute.point.quaternion(name="rotation", value=rot[type]())
            >> join
        )

    geo = GeometrySet(ensemble.object)
    assert geo.mesh

    rotation_quaternion = geo.named_attribute("rotation")
    rotation_scipy_from_gn = Rotation.from_quat(
        rotation_quaternion, scalar_first=True
    )  # blender stores quaternions as scalar-first

    # To compare the two rotation we multiply one with the inverse of the other and should get something very small
    assert (rotation_scipy * rotation_scipy_from_gn.inv()).magnitude().max() < 1e-5  # ty: ignore[unresolved-attribute]
    assert snapshot == geo


def test_load_starfiles(snapshot):
    file = data_dir / "starfile/clathrin.star"
    ensemble = mn.entities.ensemble.StarFile.load(file)
    assert ensemble._entity_type == mn.entities.base.EntityType.ENSEMBLE_STAR
    assert ensemble.props.entity_type == ensemble._entity_type.value
    assert snapshot == GeometrySet(ensemble.object)


def test_categorical_attributes(snapshot):
    file = data_dir / "starfile/cistem.star"
    ensemble = mn.entities.ensemble.StarFile.load(file)
    assert "cisTEMOriginalImageFilename" in ensemble.props.categories
    assert snapshot == GeometrySet(ensemble.object)


def test_load_ndjson_oriented(snapshot):
    file = data_dir / "cryoet/oriented_point.ndjson"
    ensemble = mn.entities.ensemble.StarFile.load(file)
    assert ensemble._entity_type == mn.entities.base.EntityType.ENSEMBLE_STAR
    assert ensemble.props.entity_type == ensemble._entity_type.value

    records = [json.loads(line) for line in open(file)]

    # positions are the voxel coordinates from the file, at world scale
    positions = databpy.named_attribute(ensemble.object, "position")
    expected = np.array(
        [[record["location"][axis] for axis in "xyz"] for record in records]
    )
    assert np.allclose(positions, expected * 0.1, atol=1e-4)

    # the stored transform combines the file's rotation matrix with the scaled
    # position, transposed on storage for Blender's column-major float4x4
    transforms = databpy.named_attribute(ensemble.object, "transform")
    expected_transforms = np.tile(np.identity(4), (len(records), 1, 1))
    expected_transforms[:, :3, :3] = [
        record["xyz_rotation_matrix"] for record in records
    ]
    expected_transforms[:, :3, 3] = expected * 0.1
    assert np.allclose(transforms, expected_transforms.transpose(0, 2, 1), atol=1e-4)

    assert snapshot == GeometrySet(ensemble.object)


def test_load_ndjson_point(snapshot):
    file = data_dir / "cryoet/point.ndjson"
    ensemble = mn.entities.ensemble.StarFile.load(file)

    records = [json.loads(line) for line in open(file)]
    positions = databpy.named_attribute(ensemble.object, "position")
    assert len(positions) == len(records)

    # plain points carry no orientation, so the stored transforms hold an
    # identity rotation with the scaled position
    transforms = databpy.named_attribute(ensemble.object, "transform")
    assert np.allclose(transforms[:, :3, :3], np.identity(3), atol=1e-4)
    assert np.allclose(transforms[:, 3, :3], positions, atol=1e-4)
    assert snapshot == GeometrySet(ensemble.object)


def test_micrograph_conversion(snapshot):
    file = data_dir / "starfile/cistem.star"
    ensemble = mn.entities.ensemble.StarFile.load(file)
    tiff_path = data_dir / "starfile/montage.tiff"
    ensemble._convert_mrc_to_tiff()
    assert tiff_path.exists()
    assert snapshot == GeometrySet(ensemble.object)
