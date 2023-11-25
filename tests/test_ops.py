import bpy
import pytest
import molecularnodes as mn
from .constants import (
    test_data_directory, 
    codes
)
from molecularnodes.io.mda import HAS_mda


if HAS_mda:
    import MDAnalysis as mda
from .utils import apply_mods, sample_attribute_to_string

# register the operators, which isn't done by default when loading bpy
# just via headless float_decimals
mn.register()


@pytest.mark.parametrize("code", codes)
def test_op_api_cartoon(snapshot, code, style = 'ribbon'):
    bpy.context.scene.MN_pdb_code = code
    bpy.context.scene.MN_import_style = style
    
    bpy.ops.mn.import_protein_rcsb()
    
    obj_1 = bpy.context.active_object
    obj_2 = mn.io.pdb.load(code, style=style)
    
    # objects being imported via each method should have identical snapshots
    for mol in [obj_1, obj_2]:
        apply_mods(mol)
        for att in mol.data.attributes.keys():
            snapshot.assert_match(
                sample_attribute_to_string(mol, att), 
                f"{att}_value.txt"
            )

def test_op_api_mda(snapshot):
    topo = str(test_data_directory / "md_ppr/box.gro")
    traj = str(test_data_directory / "md_ppr/first_5_frames.xtc")
    name = bpy.context.scene.MN_import_md_name
    
    bpy.context.scene.MN_import_md_topology  = topo
    bpy.context.scene.MN_import_md_trajectory  = traj
    
    bpy.ops.mn.import_protein_md()
    obj_1 = bpy.context.active_object
    assert obj_1.name == name
    assert not bpy.data.collections.get(f"{name}_frames")
    
    bpy.context.scene.MN_md_in_memory = True
    name = 'NewTrajectoryInMemory'
    
    obj_2, universe = mn.io.md.load(topo, traj, name = "test", style = 'ribbon', in_memory=True)
    print(list(bpy.data.objects))
    frames_coll = bpy.data.collections.get(f"{obj_2.name}_frames")
    
    assert frames_coll
    assert len(frames_coll.objects) == 5
    
    for mol in [obj_1, obj_2]:
        print(f"{mol=}")
        apply_mods(mol)
        for att in mol.data.attributes.keys():
            snapshot.assert_match(
                sample_attribute_to_string(mol, att), 
                f"{att}_value.txt"
            )
