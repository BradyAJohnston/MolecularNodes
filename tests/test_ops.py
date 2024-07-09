import bpy
import pytest
import numpy as np
import molecularnodes as mn
from molecularnodes.blender.obj import ObjectTracker, get_attribute

from .utils import sample_attribute, NumpySnapshotExtension
from .constants import data_dir, codes, attributes

# register the operators, which isn't done by default when loading bpy
# just via headless float_decimals
mn._test_register()


@pytest.mark.parametrize("code", codes)
def test_op_api_cartoon(
    snapshot_custom: NumpySnapshotExtension, code, style="ribbon", format="bcif"
):
    scene = bpy.context.scene
    scene.MN_import_node_setup = True
    scene.MN_pdb_code = code
    scene.MN_import_style = style
    scene.MN_import_node_setup = True
    scene.MN_import_build_assembly = False
    scene.MN_import_centre = False
    scene.MN_import_del_solvent = False
    scene.MN_import_format_download = format

    bpy.ops.mn.import_wwpdb()

    obj_1 = bpy.context.active_object
    obj_2 = mn.io.fetch(code, style=style, format=format, cache_dir=data_dir).object

    # objects being imported via each method should have identical snapshots
    for mol in [obj_1, obj_2]:
        for name in attributes:
            if name == "sec_struct" or name.startswith("."):
                continue
            assert snapshot_custom == sample_attribute(mol, name, evaluate=True)


@pytest.mark.parametrize("code", codes)
@pytest.mark.parametrize("file_format", ["bcif", "cif", "pdb"])
def test_op_local(snapshot_custom, code, file_format):
    scene = bpy.context.scene
    scene.MN_import_node_setup = False
    scene.MN_import_style = "spheres"
    scene.MN_import_build_assembly = False
    scene.MN_import_del_solvent = False
    scene.MN_import_format_download = file_format
    path = str(mn.io.download(code=code, format=file_format, cache=data_dir))
    scene.MN_import_local_path = path
    scene.MN_centre_type = "centroid"

    scene.MN_import_centre = False
    with ObjectTracker() as o:
        bpy.ops.mn.import_protein_local()
        bob = o.latest()

    scene.MN_import_centre = True
    with ObjectTracker() as o:
        bpy.ops.mn.import_protein_local()
        bob_centred = o.latest()

    bob_pos, bob_centred_pos = [
        sample_attribute(x, "position", evaluate=False) for x in [bob, bob_centred]
    ]

    assert snapshot_custom == bob_pos
    assert snapshot_custom == bob_centred_pos
    assert not np.allclose(bob_pos, bob_centred_pos)


def test_op_api_mda(snapshot_custom: NumpySnapshotExtension):
    topo = str(data_dir / "md_ppr/box.gro")
    traj = str(data_dir / "md_ppr/first_5_frames.xtc")
    name = bpy.context.scene.MN_import_md_name

    bpy.context.scene.MN_import_md_topology = topo
    bpy.context.scene.MN_import_md_trajectory = traj
    bpy.context.scene.MN_import_style = "ribbon"

    with ObjectTracker() as o:
        bpy.ops.mn.import_protein_md()
        obj_1 = o.latest()

    assert obj_1.name == name

    mnu = mn.io.trajectory.load(topo, traj, name="test", style="ribbon")
    obj_2 = mnu.object

    for mol in [obj_1, obj_2]:
        for att in attributes:
            assert snapshot_custom == sample_attribute(mol, att)

    # capture positions, change the frame number and test that the positions have updated
    # and cahnged
    pos_1, pos_2 = [get_attribute(x, "position") for x in [obj_1, obj_2]]
    bpy.context.scene.frame_set(4)

    assert not np.allclose(get_attribute(obj_1, "position"), pos_1)
    assert np.allclose(
        get_attribute(obj_1, "position"), get_attribute(obj_2, "position")
    )
