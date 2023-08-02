import bpy
import MolecularNodes as mn
from .utils import apply_mods, get_verts

def test_creat_obj():
    # Create a mesh object named "MyMesh" in the collection "MyCollection"
    # with vertex locations and bond edges.
    locations = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]
    bonds = [(0, 1), (1, 2), (2, 0)]
    name = "MyMesh"
    my_object = mn.obj.create_object(name, bpy.data.collections['Collection'], locations, bonds)
    
    assert len(my_object.data.vertices) == 3
    assert my_object.name == name
    assert my_object.name != "name"