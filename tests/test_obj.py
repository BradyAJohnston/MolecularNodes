import bpy
import numpy as np
import molecularnodes as mn
from molecularnodes.blender import mesh
from molecularnodes.blender import databpy as db
from .constants import data_dir

mn.register()


def test_creat_obj():
    # Create a mesh object named "MyMesh" in the collection "MyCollection"
    # with vertex locations and bond edges.
    locations = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]]
    bonds = [(0, 1), (1, 2), (2, 0)]
    name = "MyMesh"
    my_object = db.create_object(locations, bonds, name=name)

    assert len(my_object.data.vertices) == 3
    assert my_object.name == name
    assert my_object.name != "name"


def test_set_position():
    mol = mn.entities.fetch("8FAT", cache_dir=data_dir)

    pos_a = mol.named_attribute("position")

    mol.store_named_attribute(mol.named_attribute("position") + 10, name="position")

    pos_b = mol.named_attribute("position")
    print(f"{pos_a=}")
    print(f"{pos_b=}")

    assert not np.isclose(pos_a, pos_b).all()
    assert np.isclose(pos_a, pos_b - 10, rtol=0.1).all()


def test_eval_mesh():
    a = db.create_object(np.zeros((3, 3)))
    assert len(a.data.vertices) == 3
    b = db.create_object(np.zeros((5, 3)))
    assert len(b.data.vertices) == 5
    assert len(mesh.evaluate_using_mesh(b).data.vertices) == 5


def test_matrix_read_write():
    obj = db.create_object(np.zeros((5, 3)))
    arr = np.array((5, 4, 4), float)
    arr = np.random.rand(5, 4, 4)

    db.store_named_attribute(
        obj=obj, data=arr, name="test_matrix", atype=db.AttributeTypes.FLOAT4X4
    )

    assert np.allclose(mesh.named_attribute(obj, "test_matrix"), arr)

    arr2 = np.random.rand(5, 4, 4)
    db.store_named_attribute(obj=obj, data=arr2, name="test_matrix2")

    assert not np.allclose(mesh.named_attribute(obj, "test_matrix2"), arr)
