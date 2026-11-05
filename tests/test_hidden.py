"""Atoms with the `hidden` attribute are left out of every style."""

import bpy
import numpy as np
import pytest
import molecularnodes as mn
from molecularnodes.blender import utils as blender_utils
from molecularnodes.nodes import geometry as mg
from .constants import data_dir

STYLES = ["spheres", "sticks", "ball_and_stick", "cartoon", "ribbon", "surface"]


def _molecule():
    return mn.Molecule.fetch("4ozs", cache=data_dir, format="bcif")


def _hidden_mask(mol):
    res_id = mol.named_attribute("res_id")
    return res_id < np.median(res_id)


def evaluated_vertex_count(obj):
    depsgraph = bpy.context.evaluated_depsgraph_get()
    return len(obj.evaluated_get(depsgraph).data.vertices)


def _rendered(mol):
    points = blender_utils.evaluated_points(mol.object)
    return points[np.lexsort(points.T)]


@pytest.fixture(autouse=True)
def canvas():
    return mn.Canvas(mn.scene.Cycles(samples=1, device="CPU"))


@pytest.mark.parametrize("style", STYLES)
def test_hidden_matches_selection(style):
    hidden = _molecule()
    mask = _hidden_mask(hidden)
    assert 0 < mask.sum() < len(mask)
    hidden.store_named_attribute(mask, "hidden", atype="BOOLEAN")
    hidden.add_style(style)

    selected = _molecule()
    selected.store_named_attribute(~mask, "shown", atype="BOOLEAN")
    selected.add_style(style, selection="shown")

    assert np.allclose(_rendered(hidden), _rendered(selected))


@pytest.mark.parametrize("style", STYLES)
def test_include_hidden(style):
    plain = _molecule()
    plain.add_style(style)

    hidden = _molecule()
    hidden.store_named_attribute(_hidden_mask(hidden), "hidden", atype="BOOLEAN")
    hidden.add_style(style, include_hidden=True)

    assert np.allclose(_rendered(hidden), _rendered(plain))


def test_hidden_updates_existing_styles():
    mol = _molecule()
    mol.add_style("spheres")
    before = len(_rendered(mol))
    mol.store_named_attribute(_hidden_mask(mol), "hidden", atype="BOOLEAN")
    assert len(_rendered(mol)) < before


def test_hidden_wins_over_selection():
    hidden = _molecule()
    mask = _hidden_mask(hidden)
    hidden.store_named_attribute(mask, "hidden", atype="BOOLEAN")
    hidden.store_named_attribute(np.ones_like(mask), "wanted", atype="BOOLEAN")
    hidden.add_style("spheres", selection="wanted")

    selected = _molecule()
    selected.store_named_attribute(~mask, "shown", atype="BOOLEAN")
    selected.add_style("spheres", selection="shown")

    assert np.allclose(_rendered(hidden), _rendered(selected))


def test_hidden_atoms_are_not_computed_by_style_hbonds(fetch):
    mol = fetch("2M6Q")
    obj = mol.calculations.hbonds().object
    n_bonds = len(obj.data.edges)
    assert n_bonds > 0

    def vertices(**inputs):
        with mol.tree.reset() as (atoms, join):
            atoms >> mg.StyleHBonds(data_object=obj, quality=2, **inputs) >> join
        return evaluated_vertex_count(mol.object)

    visible = vertices()
    assert visible > 0
    hydrogens = mol.named_attribute("atomic_number") == 1
    mol.store_named_attribute(hydrogens, "hidden", atype="BOOLEAN")
    assert vertices() == 0
    assert vertices(include_hidden=True) == visible


def test_hidden_in_oxdna_style():
    mol = mn.entities.OXDNA.load(
        data_dir / "oxdna/linear.top", data_dir / "oxdna/linear_traj.dat"
    )
    before = len(_rendered(mol))
    mask = np.arange(len(mol.position)) < len(mol.position) // 2
    mol.store_named_attribute(mask, "hidden", atype="BOOLEAN")
    assert len(_rendered(mol)) < before
