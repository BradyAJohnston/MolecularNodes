import bpy
import MolecularNodes as mn

from .utils import get_verts

# register the operators, which isn't done by default when loading bpy
# just via headless python
mn.register()

def compare_op_api(code, style = "atoms", apply = True):
    bpy.context.scene.MN_pdb_code = code
    bpy.context.scene.MN_import_default_style = style
    
    bpy.ops.mn.import_protein_rcsb()
    obj_1 = bpy.data.objects[code]
    obj_2 = mn.load.molecule_rcsb(code, starting_style=style)
    
    v1 = get_verts(obj_1, apply_modifiers=apply)
    v2 = get_verts(obj_2, apply_modifiers=apply)
    return  v1 == v2

def test_op_rcsb_6n2y():
    assert compare_op_api('6n2y', style = "cartoon", apply = True)

def test_op_rcsb_4ozs():
    assert compare_op_api('4ozs', style = "atoms", apply = False)