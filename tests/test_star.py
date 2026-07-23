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
    ensemble = mn.entities.ensemble.load_starfile(file)

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
    ensemble = mn.entities.ensemble.load_starfile(file)
    assert ensemble._entity_type == mn.entities.base.EntityType.ENSEMBLE_STAR
    assert ensemble.object.mn.entity_type == ensemble._entity_type.value
    assert snapshot == GeometrySet(ensemble.object)


def test_categorical_attributes(snapshot):
    file = data_dir / "starfile/cistem.star"
    ensemble = mn.entities.ensemble.load_starfile(file)
    assert "cisTEMOriginalImageFilename_categories" in ensemble.object
    assert snapshot == GeometrySet(ensemble.object)


def test_micrograph_conversion(snapshot):
    file = data_dir / "starfile/cistem.star"
    ensemble = mn.entities.ensemble.load_starfile(file)
    tiff_path = data_dir / "starfile/montage.tiff"
    ensemble._convert_mrc_to_tiff()
    assert tiff_path.exists()
    assert snapshot == GeometrySet(ensemble.object)
