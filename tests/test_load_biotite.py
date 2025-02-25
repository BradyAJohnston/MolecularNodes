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

    # use the class method
    mol = mn.entities.Molecule.from_array(arr)
    assert isinstance(mol,  mn.entities.Molecule)
    assert mol.file_path == None
    assert mol.file == "ARRAY_LOADED_DIRECTLY"
    assert mol._frames_collection == None
    #assert mol._entity_type == EntityType.MOLECULE
    assert mol._assemblies() == None
    assert mol.n_atoms == 4730
