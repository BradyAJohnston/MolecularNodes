import bpy
import pytest
import starfile
import importlib
from scipy.spatial.transform import Rotation as R

import molecularnodes as mn

from .constants import data_dir

mn.unregister()
mn.register()


@pytest.mark.parametrize("type", ["cistem", "relion"])
def test_starfile_attributes(type):
    file = data_dir / f"{type}.star"
    ensemble = mn.io.star.load(file)

    star = starfile.read(file)

    if type == "relion":
        df = star["particles"].merge(star["optics"], on="rlnOpticsGroup")
        euler_angles = df[["rlnAngleRot", "rlnAngleTilt", "rlnAnglePsi"]].to_numpy()

    elif type == "cistem":
        df = star
        euler_angles = df[["cisTEMAnglePhi", "cisTEMAngleTheta", "cisTEMAnglePsi"]].to_numpy()

    # Calculate Scipy rotation from the euler angles
    rot_from_euler = R.from_euler(seq="ZYZ", angles=euler_angles, degrees=True).inv()

    # Activate the rotation debug mode in the nodetreee and get the quaternion attribute
    debugnode = mn.blender.nodes.star_node(ensemble.node_group).node_tree.nodes["Switch.001"]
    debugnode.inputs["Switch"].default_value = True
    quat_attribute = ensemble.get_attribute("MNDEBUGEuler", evaluate=True)

    # Convert from blender to scipy conventions and then into Scipy rotation
    rot_from_geo_nodes = R.from_quat(quat_attribute[:, [1, 2, 3, 0]])

    # To compare the two rotation with multiply one with the inverse of the other
    assert (rot_from_euler * rot_from_geo_nodes.inv()).magnitude().max() < 1e-5


def test_categorical_attributes():
    file = data_dir / "cistem.star"
    ensemble = mn.io.star.load(file)
    assert "cisTEMOriginalImageFilename_categories" in ensemble.object


def test_micrograph_conversion():
    file = data_dir / "cistem.star"
    ensemble = mn.io.star.load(file)
    tiff_path = data_dir / "montage.tiff"
    tiff_path.unlink(missing_ok=True)
    ensemble._convert_mrc_to_tiff()
    assert tiff_path.exists()


def test_micrograph_loading():
    import bpy

    file = data_dir / "cistem.star"
    tiff_path = data_dir / "montage.tiff"
    tiff_path.unlink(missing_ok=True)

    ensemble = mn.io.star.load(file)
    assert not tiff_path.exists()
    ensemble.star_node.inputs["Show Micrograph"].default_value = True
    bpy.context.evaluated_depsgraph_get().update()
    assert tiff_path.exists()
    # Ensure montage get only loaded once
    assert sum(1 for image in bpy.data.images.keys() if "montage" in image) == 1
    assert ensemble.micrograph_material.node_tree.nodes["Image Texture"].image.name == "montage.tiff"
    assert ensemble.star_node.inputs["Micrograph"].default_value.name == "montage.tiff"


@pytest.mark.skipif(importlib.util.find_spec("pyopenvdb"), reason="Test may segfault on GHA")
def test_rehydration(tmp_path):
    bpy.ops.wm.read_homefile()
    ensemble = mn.io.star.load(data_dir / "cistem.star")
    bpy.ops.wm.save_as_mainfile(filepath=str(tmp_path / "test.blend"))
    assert ensemble._update_micrograph_texture in bpy.app.handlers.depsgraph_update_post
    bpy.ops.wm.read_homefile()
    assert ensemble._update_micrograph_texture not in bpy.app.handlers.depsgraph_update_post
    bpy.ops.wm.open_mainfile(filepath=str(tmp_path / "test.blend"))
    new_ensemble = bpy.types.Scene.MN_starfile_ensembles[0]
    assert new_ensemble._update_micrograph_texture in bpy.app.handlers.depsgraph_update_post
    assert new_ensemble.data.equals(ensemble.data)
