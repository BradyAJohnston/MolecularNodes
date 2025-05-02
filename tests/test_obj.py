import databpy
import numpy as np
import molecularnodes as mn
from molecularnodes.blender import mesh
from .constants import data_dir


def test_set_position():
    mol = mn.Molecule.fetch("8FAT", cache=data_dir)
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
