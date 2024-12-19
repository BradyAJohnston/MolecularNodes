import numpy as np
import molecularnodes as mn
from molecularnodes.blender import mesh
import databpy
from databpy import BlenderObject
from databpy.object import LinkedObjectError
from .constants import data_dir
import bpy
import pytest


def test_creat_obj():
    # Create a mesh object named "MyMesh" in the collection "MyCollection"
    # with vertex locations and bond edges.
    locations = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]
    bonds = [(0, 1), (1, 2), (2, 0)]
    name = "MyMesh"
    my_object = databpy.create_object(locations, bonds, name=name)

    assert len(my_object.data.vertices) == 3
    assert my_object.name == name
    assert my_object.name != "name"


def test_BlenderObject():
    bob = BlenderObject(None)

    with pytest.raises(LinkedObjectError):
        bob.object
    with pytest.raises(LinkedObjectError):
        bob.name
    with pytest.raises(LinkedObjectError):
        bob.name = "testing"

    bob = BlenderObject(bpy.data.objects["Cube"])
    assert bob.name == "Cube"
    bob.name = "NewName"
    with pytest.raises(KeyError):
        bpy.data.objects["Cube"]
    assert bob.name == "NewName"


def test_bob():
    mol = mn.entities.fetch("8H1B", cache_dir=data_dir)
    assert isinstance(mol, BlenderObject)
    with pytest.raises(NotImplementedError):
        mol.set_frame(10)

    with pytest.raises(ValueError):
        mol.frames

    mol2 = mn.entities.fetch("1NMR", cache_dir=data_dir)
    assert isinstance(mol2.frames, bpy.types.Collection)
    assert mol2.name == "1NMR"


def test_set_position():
    mol = mn.entities.fetch("8FAT", cache_dir=data_dir)
    pos_a = mol.position
    mol.position += 10
    pos_b = mol.position
    assert not np.allclose(pos_a, pos_b)
    assert np.allclose(pos_a, pos_b - 10, rtol=0.1)


def test_change_names():
    bob_cube = databpy.BlenderObject("Cube")
    assert bob_cube.name == "Cube"
    with databpy.ObjectTracker() as o:
        bpy.ops.mesh.primitive_cylinder_add()
        bob_cyl = databpy.BlenderObject(o.latest())

    assert bob_cyl.name == "Cylinder"
    assert len(bob_cube) != len(bob_cyl)

    # rename the objects, but separately to the linked BlenderObject, so that the
    # reference will have to be rebuilt from the .uuid when the names don't match
    bpy.data.objects["Cylinder"].name = "Cylinder2"
    bpy.data.objects["Cube"].name = "Cylinder"

    # ensure that the reference to the actul object is updated, so that even if the name has
    # changed the reference is reconnected via the .uuid
    assert len(bob_cube) == 8
    assert bob_cube.name == "Cylinder"
    assert bob_cyl.name == "Cylinder2"


def test_eval_mesh():
    a = databpy.create_object(np.zeros((3, 3)))
    assert len(a.data.vertices) == 3
    b = databpy.create_object(np.zeros((5, 3)))
    assert len(b.data.vertices) == 5
    assert len(mesh.evaluate_using_mesh(b).data.vertices) == 5


def test_matrix_read_write():
    bob = databpy.create_bob(np.zeros((5, 3)))
    arr = np.array((5, 4, 4), float)
    arr = np.random.rand(5, 4, 4)

    bob.store_named_attribute(
        data=arr, name="test_matrix", atype=databpy.AttributeTypes.FLOAT4X4
    )

    assert np.allclose(bob.named_attribute("test_matrix"), arr)
    arr2 = np.random.rand(5, 4, 4)
    bob.store_named_attribute(data=arr2, name="test_matrix2")
    assert (
        bob.object.data.attributes["test_matrix2"].data_type
        == databpy.AttributeTypes.FLOAT4X4.value.type_name
    )
    assert not np.allclose(bob.named_attribute("test_matrix2"), arr)
