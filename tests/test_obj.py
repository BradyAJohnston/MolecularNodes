import numpy as np
import molecularnodes as mn
from molecularnodes.blender import mesh
import databpy
from databpy import BlenderObject
from databpy.object import LinkedObjectError
from .constants import data_dir
import bpy
import pytest


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


def test_eval_mesh():
    a = databpy.create_object(np.zeros((3, 3)))
    assert len(a.data.vertices) == 3
    b = databpy.create_object(np.zeros((5, 3)))
    assert len(b.data.vertices) == 5
    assert len(mesh.evaluate_using_mesh(b).data.vertices) == 5
