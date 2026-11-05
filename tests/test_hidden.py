"""Atoms with the `hidden` attribute are left out of every style."""

import warnings
import bpy
import MDAnalysis as mda
import numpy as np
import pytest
from databpy import ObjectTracker
import molecularnodes as mn
from molecularnodes.blender import utils as blender_utils
from molecularnodes.nodes import geometry as mg
from molecularnodes.session import get_session
from .constants import data_dir
from .utils import mn_ops

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


def test_hidden_is_created_on_first_read():
    mol = _molecule()
    assert "hidden" not in mol.list_attributes()
    assert not mol.hidden.any()
    assert mol.hidden.dtype == bool
    assert "hidden" in mol.list_attributes()


def test_hidden_writes_back():
    mol = _molecule()
    mol.hidden[:10] = True
    assert mol.named_attribute("hidden").sum() == 10
    assert len(mol.universe.atoms[np.asarray(mol.hidden)]) == 10

    mol.hidden = np.arange(len(mol.position)) < 3
    assert mol.named_attribute("hidden").sum() == 3


def test_hide_adds_and_show_removes():
    mol = _molecule()
    ca = mol.universe.select_atoms("name CA")
    assert mol.hide("name CA") is mol
    assert mol.hidden.sum() == len(ca)
    mol.hide("name N")
    n_both = len(mol.universe.select_atoms("name CA N"))
    assert mol.hidden.sum() == n_both
    assert mol.show("name N") is mol
    assert mol.hidden.sum() == len(ca)
    mol.show()
    assert not mol.hidden.any()


def test_hide_requires_a_selection():
    with pytest.raises(TypeError):
        _molecule().hide()


def test_hide_accepts_each_selection_type():
    mol = _molecule()
    expected = mol.universe.select_atoms("name CA")
    mask = np.zeros(len(mol.position), dtype=bool)
    mask[expected.ix] = True
    mol.store_named_attribute(mask, "alpha", atype="BOOLEAN")
    for selection in ("name CA", "alpha", expected, mask):
        mol.show().hide(selection)
        assert np.array_equal(mol.hidden, mask)


def test_hide_rejects_what_cannot_be_evaluated_once():
    mol = _molecule()
    updating = mol.universe.select_atoms("around 5 name CA", updating=True)
    with pytest.raises(ValueError, match="custom selection"):
        mol.hide(updating)
    with pytest.raises(ValueError, match="custom selection"):
        mol.show(updating)
    with pytest.raises(TypeError, match="callable"):
        mol.hide(lambda: mg.IsHydrogen())
    with pytest.raises(ValueError, match="Invalid selection"):
        mol.hide("not a ( selection")
    with pytest.raises(ValueError, match="one value for each"):
        mol.hide(np.array([True, False]))


@pytest.mark.parametrize("method", ["from_string", "from_atomgroup"])
def test_hidden_is_a_reserved_selection_name(method):
    mol = _molecule()
    selection = "name CA" if method == "from_string" else mol.universe.atoms[:5]
    with pytest.raises(ValueError, match="reserved"):
        getattr(mol.selections, method)(selection, name="hidden")


def test_load_and_fetch_hide():
    path = data_dir / "1ATN.pdb"
    mol = mn.Molecule.load(path, hide="name CA", style="spheres")
    assert mol.hidden.sum() == len(mol.universe.select_atoms("name CA"))
    assert mol.props.hidden_selection == "name CA"

    fetched = mn.Molecule.fetch("4ozs", cache=data_dir, format="bcif", hide="name CA")
    assert fetched.hidden.sum() == len(fetched.universe.select_atoms("name CA"))


def test_hidden_selection_text_replaces_hidden_atoms():
    mol = _molecule()
    mol.hide("name N")
    mol.props.hidden_selection = "name CA"
    assert mol.hidden.sum() == len(mol.universe.select_atoms("name CA"))

    mol.props.hidden_selection = "not a ( selection"
    assert "Invalid selection" in mol.props.hidden_message
    assert mol.hidden.sum() == len(mol.universe.select_atoms("name CA"))

    mol.props.hidden_selection = ""
    assert mol.props.hidden_message == ""
    assert not mol.hidden.any()


def test_python_edits_clear_the_selection_text():
    mol = _molecule()
    mol.props.hidden_selection = "name CA"
    mol.hide("name N")
    assert mol.props.hidden_selection == ""
    n_hidden = len(mol.universe.select_atoms("name CA N"))
    assert mol.hidden.sum() == n_hidden


def test_add_style_warns_when_selection_is_all_hidden():
    mol = _molecule()
    mol.hide("name CA")
    with pytest.warns(UserWarning, match="is hidden"):
        mol.add_style("spheres", selection="name CA")
    with warnings.catch_warnings():
        warnings.filterwarnings("error", message=".*is hidden.*")
        mol.add_style("spheres", selection="name CA", include_hidden=True)
        mol.add_style("spheres", selection="name CA or name N")


def test_calculations_still_see_hidden_atoms():
    mol = _molecule()
    n_atoms = len(mol.universe.atoms)
    mol.hide("all")
    assert len(mol.universe.select_atoms("all")) == n_atoms
    assert len(mol.position) == n_atoms


def test_import_operator_hide():
    with ObjectTracker() as tracker:
        mn_ops.import_molecule(
            code="4ozs",
            file_format="bcif",
            cache_dir=str(data_dir),
            style="spheres",
            hide="name CA",
        )
    mol = get_session().match(tracker.latest())
    assert isinstance(mol, mn.Molecule)
    assert mol.hidden.sum() == len(mol.universe.select_atoms("name CA"))
    assert mol.props.hidden_selection == "name CA"


def test_import_operator_hides_water_and_hydrogens_by_default():
    with ObjectTracker() as tracker:
        mn_ops.import_molecule(
            code="4ozs", file_format="bcif", cache_dir=str(data_dir), style="spheres"
        )
    mol = get_session().match(tracker.latest())
    assert isinstance(mol, mn.Molecule)
    assert mol.hidden.sum() == len(mol.universe.select_atoms("water or element H")) > 0


def test_import_operator_default_without_elements():
    with ObjectTracker() as tracker:
        mn_ops.import_molecule(
            method="local",
            filepath=str(data_dir / "md_ppr/box.gro"),
            trajectory=str(data_dir / "md_ppr/first_5_frames.xtc"),
        )
    mol = get_session().match(tracker.latest())
    assert isinstance(mol, mn.Molecule)
    assert mol.props.hidden_selection == "water"


def test_import_operator_can_hide_nothing():
    with ObjectTracker() as tracker:
        mn_ops.import_molecule(
            code="4ozs",
            file_format="bcif",
            cache_dir=str(data_dir),
            style="spheres",
            hide="",
        )
    mol = get_session().match(tracker.latest())
    assert isinstance(mol, mn.Molecule)
    assert not mol.hidden.any()


def test_hidden_survives_frame_changes():
    universe = mda.Universe(
        data_dir / "md_ppr/box.gro", data_dir / "md_ppr/first_5_frames.xtc"
    )
    mol = mn.Molecule(universe)
    mol.add_style("spheres")
    mol.hide("index 0:99")
    expected = np.array(mol.hidden)
    assert expected.any()
    mol.set_frame(3)
    assert np.array_equal(mol.hidden, expected)
