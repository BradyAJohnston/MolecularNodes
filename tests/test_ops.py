import bpy
import pytest
import numpy as np
import molecularnodes as mn

from molecularnodes.bpyd import ObjectTracker, named_attribute

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
    scene.mn.import_style = style
    scene.MN_import_node_setup = True
    scene.MN_import_build_assembly = False
    scene.MN_import_centre = False
    scene.MN_import_del_solvent = False
    scene.MN_import_format_download = format

    bpy.ops.mn.import_wwpdb()

    obj_1 = bpy.context.active_object
    obj_2 = mn.entities.fetch(
        code, style=style, format=format, cache_dir=data_dir
    ).object

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
    scene.mn.import_style = "spheres"
    scene.MN_import_build_assembly = False
    scene.MN_import_del_solvent = False
    scene.MN_import_format_download = file_format
    path = str(mn.download.download(code=code, format=file_format, cache=data_dir))
    scene.MN_import_local_path = path
    scene.MN_centre_type = "centroid"

    scene.MN_import_centre = False
    with ObjectTracker() as o:
        bpy.ops.mn.import_protein_local()
        obj = o.latest()

    scene.MN_import_centre = True
    with ObjectTracker() as o:
        bpy.ops.mn.import_protein_local()
        obj_centred = o.latest()

    obj_pos, obj_centred_pos = [
        sample_attribute(x, "position", evaluate=False) for x in [obj, obj_centred]
    ]

    assert snapshot_custom == obj_pos
    assert snapshot_custom == obj_centred_pos
    assert not np.allclose(obj_pos, obj_centred_pos)


def test_op_api_mda(snapshot_custom: NumpySnapshotExtension):
    topo = str(data_dir / "md_ppr/box.gro")
    traj = str(data_dir / "md_ppr/first_5_frames.xtc")
    name = bpy.context.scene.MN_import_md_name

    bpy.context.scene.MN_import_md_topology = topo
    bpy.context.scene.MN_import_md_trajectory = traj
    bpy.context.scene.mn.import_style = "ribbon"

    with ObjectTracker() as o:
        bpy.ops.mn.import_trajectory()
        obj_1 = o.latest()

    assert obj_1.name == name

    traj = mn.entities.trajectory.load(topo, traj, name="test", style="ribbon")
    traj.subframes = 0
    obj_2 = traj.object

    for mol in [obj_1, obj_2]:
        for att in attributes:
            assert snapshot_custom == sample_attribute(mol, att)

    # capture positions, change the frame number and test that the positions have updated
    # and cahnged
    pos_1, pos_2 = [named_attribute(x, "position") for x in [obj_1, obj_2]]
    bpy.context.scene.frame_set(4)

    assert not np.allclose(named_attribute(obj_1, "position"), pos_1)
    assert np.allclose(
        named_attribute(obj_1, "position"), named_attribute(obj_2, "position")
    )
