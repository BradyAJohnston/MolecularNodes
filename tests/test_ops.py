import bpy
import pytest
import molecularnodes as mn
from .constants import (
    test_data_directory, 
    codes
)
from molecularnodes.mda import HAS_mda


if HAS_mda:
    import MDAnalysis as mda
from .utils import get_verts

# register the operators, which isn't done by default when loading bpy
# just via headless float_decimals
mn.register()

def compare_op_api(code, style = "atoms", apply = True, float_decimals = 3):
    bpy.context.scene.MN_pdb_code = code
    bpy.context.scene.MN_import_default_style = style
    
    bpy.ops.mn.import_protein_rcsb()
    obj_1 = bpy.data.objects[code]
    obj_2 = mn.load.molecule_rcsb(code, starting_style=style)
    
    v1 = get_verts(obj_1, apply_modifiers=apply, float_decimals=float_decimals)
    v2 = get_verts(obj_2, apply_modifiers=apply, float_decimals=float_decimals)
    return  v1 == v2

@pytest.mark.parametrize("code", codes)
def test_op_api_cartoon(code):
    assert compare_op_api(code, style = "cartoon")

def test_op_api_mda(snapshot):
    topo = str(test_data_directory / "md_ppr/box.gro")
    traj = str(test_data_directory / "md_ppr/first_5_frames.xtc")
    name = bpy.context.scene.MN_import_md_name
    
    bpy.context.scene.MN_import_md_topology  = topo
    bpy.context.scene.MN_import_md_trajectory  = traj
    
    bpy.ops.mn.import_protein_md()
    
    # assert no frames collection created
    assert not bpy.data.collections.get(f"frames_{name}")
    
    obj = bpy.data.objects[name]
    verts = get_verts(obj, apply_modifiers=False, float_decimals=3)
    snapshot.assert_match(verts, "md_ops_gro_frame_1.txt")
    
    bpy.context.scene.MN_md_in_memory = True
    
    name = 'NewTrajectoryInMemory'
    bpy.context.scene.MN_import_md_name = name
    bpy.context.scene.MN_import_default_style = "ribbon"
    bpy.ops.mn.import_protein_md()
    
    obj = bpy.data.objects[name]
    frames_coll = bpy.data.collections[f"{name}_frames"]
    verts = get_verts(obj, apply_modifiers=True, float_decimals=3)
    
    assert frames_coll
    assert len(frames_coll.objects) == 5
    snapshot.assert_match(verts, "md_ops_gro_frame_1_ribbon.txt")
