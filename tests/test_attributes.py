import itertools
import numpy as np
import pytest
import molecularnodes as mn
from molecularnodes.nodes.geometry import StyleRibbon, StyleSpheres
from .constants import codes, data_dir
from .utils import GeometrySet

formats = ["pdb", "cif", "bcif"]


@pytest.mark.parametrize("code, format", itertools.product(codes, formats))
def test_attribute(snapshot, code, format):
    mol = mn.Molecule.fetch(code, cache=data_dir, format=format)
    assert snapshot == GeometrySet(mol.object)


def test_store_named_attribute(snapshot_custom):
    mol = mn.Molecule.fetch("8H1B", cache=data_dir, format="bcif")
    before = mol.named_attribute("position")
    mol.store_named_attribute(mol.named_attribute("position") + 10, "position")
    after = mol.named_attribute("position")

    assert not np.allclose(before, after)


def test_uv_map(snapshot_custom):
    mol = mn.Molecule.fetch("1cd3", cache=data_dir, format="bcif")
    with mol.tree.reset() as (atoms, join):
        atoms >> StyleRibbon(uv_map=True, quality=1) >> join
    assert snapshot_custom == mol.named_attribute("uv_map", evaluate=True)[:1000]
    assert snapshot_custom == mol.named_attribute("uv_map", evaluate=True)[-1000:]


def test_bond_attributes(snapshot):
    mol = mn.Molecule.fetch("1BNA", cache=data_dir, format="bcif")
    with mol.tree.reset() as (atoms, join):
        atoms >> StyleSpheres(geometry="Mesh") >> join

    assert snapshot == GeometrySet(mol.object)
