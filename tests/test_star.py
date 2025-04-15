import pytest
import starfile
from scipy.spatial.transform import Rotation
import molecularnodes as mn
from .constants import data_dir


@pytest.mark.parametrize("type", ["cistem", "relion"])
def test_starfile_attributes(type):
    file = data_dir / f"starfile/{type}.star"
    ensemble = mn.entities.ensemble.load_starfile(file)

    star = starfile.read(file)

    if type == "relion":
        df = star["particles"].merge(star["optics"], on="rlnOpticsGroup")
        euler_angles = df[["rlnAngleRot", "rlnAngleTilt", "rlnAnglePsi"]].to_numpy()

    elif type == "cistem":
        df = star
        euler_angles = df[
            ["cisTEMAnglePhi", "cisTEMAngleTheta", "cisTEMAnglePsi"]
        ].to_numpy()

    # Calculate Scipy rotation from the euler angles
    # Note: rot_from_euler = quats
    rot_from_euler = Rotation.from_euler(
        seq="ZYZ", angles=euler_angles, degrees=True
    ).inv()

    # Convert from blender to scipy conventions and then into Scipy rotation
    quat_attribute = ensemble.named_attribute("rotation")
    rot_from_geo_nodes = Rotation.from_quat(quat_attribute[:, [1, 2, 3, 0]])

    # To compare the two rotation with multiply one with the inverse of the other
    assert (rot_from_euler * rot_from_geo_nodes.inv()).magnitude().max() < 1e-5


def test_load_starfiles():
    file = data_dir / "starfile/clathrin.star"
    _ensemble = mn.entities.ensemble.load_starfile(file)


def test_categorical_attributes():
    file = data_dir / "starfile/cistem.star"
    ensemble = mn.entities.ensemble.load_starfile(file)
    assert "cisTEMOriginalImageFilename_categories" in ensemble.object


def test_micrograph_conversion():
    file = data_dir / "starfile/cistem.star"
    ensemble = mn.entities.ensemble.load_starfile(file)
    tiff_path = data_dir / "starfile/montage.tiff"
    ensemble._convert_mrc_to_tiff()
    assert tiff_path.exists()


def test_micrograph_loading():
    import bpy

    file = data_dir / "starfile/cistem.star"
    tiff_path = data_dir / "starfile/montage.tiff"
    ensemble = mn.entities.ensemble.load_starfile(file)
    ensemble.star_node.inputs["Show Micrograph"].default_value = True
    bpy.context.evaluated_depsgraph_get().update()
    assert tiff_path.exists()
    # Ensure montage get only loaded once
    assert sum(1 for image in bpy.data.images.keys() if "montage" in image) == 1
    assert (
        ensemble.micrograph_material.node_tree.nodes["Image Texture"].image.name
        == "montage.tiff"
    )
    assert ensemble.star_node.inputs["Micrograph"].default_value.name == "montage.tiff"
