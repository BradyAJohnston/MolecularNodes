import bpy
import pytest
import molecularnodes as mn

from .utils import sample_attribute, NumpySnapshotExtension
from .constants import (
    data_dir,
    codes,
    attributes
)

# register the operators, which isn't done by default when loading bpy
# just via headless float_decimals
mn.unregister()
mn.register()


@pytest.mark.parametrize("code", codes)
def test_op_api_cartoon(snapshot: NumpySnapshotExtension, code, style='ribbon', format="bcif"):
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
    obj_2 = mn.io.fetch(code, style=style, format=format).object

    # objects being imported via each method should have identical snapshots
    for mol in [obj_1, obj_2]:
        for name in attributes:
            if name == "sec_struct" or name.startswith("."):
                continue
            assert snapshot == sample_attribute(mol, name, evaluate=True)


def test_op_api_mda(snapshot: NumpySnapshotExtension):
    topo = str(data_dir / "md_ppr/box.gro")
    traj = str(data_dir / "md_ppr/first_5_frames.xtc")
    name = bpy.context.scene.MN_import_md_name

    bpy.context.scene.MN_import_md_topology = topo
    bpy.context.scene.MN_import_md_trajectory = traj
    bpy.context.scene.MN_import_style = 'ribbon'

    bpy.ops.mn.import_protein_md()
    obj_1 = bpy.context.active_object
    assert obj_1.name == name
    assert not bpy.data.collections.get(f"{name}_frames")

    bpy.context.scene.MN_md_in_memory = True
    name = 'NewTrajectoryInMemory'

    obj_2, universe = mn.io.md.load(topo, traj, name="test", style='ribbon')
    frames_coll = bpy.data.collections.get(f"{obj_2.name}_frames")

    assert not frames_coll

    for mol in [obj_1, obj_2]:
        for att in attributes:
            assert snapshot == sample_attribute(mol, att)
