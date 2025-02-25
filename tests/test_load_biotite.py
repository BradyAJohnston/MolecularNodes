import bpy
import numpy as np
import pytest
import itertools
import molecularnodes as mn
import databpy as db
from biotite.structure import io
from biotite.structure import AtomArray, AtomArrayStack
from .constants import data_dir, codes, attributes
from .utils import NumpySnapshotExtension




def test_loading():
    test_file =  data_dir / "1f2n.bcif"
    arr = io.load_structure(test_file, template=None)

    assert isinstance(arr, AtomArray)

    # method one. overrid the __init__
    # mol = mn.entities.Molecule(arr)
    # assert isinstance(mol,  mn.entities.Molecule)

    # use the class method
    mol = mn.entities.Molecule.from_array(arr)
    assert isinstance(mol,  mn.entities.Molecule)

    mol = mn.entities.Molecule.load(arr)
    assert isinstance(mol,  mn.entities.Molecule)
