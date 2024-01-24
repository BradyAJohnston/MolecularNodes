import molecularnodes as mn
import pytest
import numpy as np
from scipy.spatial.transform import Rotation as R
import starfile
from .constants import data_dir

mn.unregister()
mn.register()


@pytest.mark.parametrize("type", ["cistem", "relion"])
def test_starfile_attributes(type, snapshot):
    file = data_dir / f"{type}.star"
    obj = mn.io.star.load(file)

    star = starfile.read(file)

    if type == 'relion':
        df = star['particles'].merge(star['optics'], on='rlnOpticsGroup')
        euler_angles = df[['rlnAngleRot',
                           'rlnAngleTilt', 'rlnAnglePsi']].to_numpy()

    elif type == 'cistem':
        df = star
        euler_angles = df[['cisTEMAnglePhi',
                           'cisTEMAngleTheta', 'cisTEMAnglePsi']].to_numpy()

    # Calculate Scipy rotation from the euler angles
    rot_from_euler = quats = R.from_euler(
        seq='ZYZ', angles=euler_angles, degrees=True
    ).inv()

    # Activate the rotation debug mode in the nodetreee and get the quaternion attribute
    star_node_group = mn.blender.nodes.get_star_node(obj)
    debugnode = star_node_group.node_tree.nodes['Switch.001']
    debugnode.inputs[1].default_value = True
    evaluated = obj.evaluated_get(mn.bpy.context.evaluated_depsgraph_get())
    quat_attribute = mn.blender.obj.get_attribute(evaluated, 'MNDEBUGEuler')

    # Convert from blender to scipy conventions and then into Scipy rotation
    rot_from_geo_nodes = R.from_quat(quat_attribute[:, [1, 2, 3, 0]])

    # To compare the two rotation with multiply one with the inverse of the other
    assert (rot_from_euler * rot_from_geo_nodes.inv()).magnitude().max() < 1e-5
