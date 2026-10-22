import pickle
import bpy
import MDAnalysis as mda
import numpy as np
import pytest
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from nodebpy import geometry as g
import molecularnodes as mn
from molecularnodes.nodes.geometry import StyleHBonds
from .constants import data_dir

pytestmark = [
    pytest.mark.filterwarnings("ignore:.*Empty string to select atoms.*:UserWarning"),
    pytest.mark.filterwarnings(
        "ignore:.*Unknown masses are set to.*:PendingDeprecationWarning"
    ),
]


@pytest.fixture
def universe():
    # gro topology: no bonds or elements, so the name and distance based
    # defaults are used
    return mda.Universe(
        data_dir / "md_ppr/box.gro", data_dir / "md_ppr/first_5_frames.xtc"
    )


def edges(obj: bpy.types.Object) -> np.ndarray:
    pairs = np.zeros(len(obj.data.edges) * 2, dtype=int)
    obj.data.edges.foreach_get("vertices", pairs)
    return pairs.reshape(-1, 2)


def expected_pairs(universe: mda.Universe, frame: int, **kwargs) -> np.ndarray:
    hba = HydrogenBondAnalysis(universe, update_selections=False, **kwargs)
    hba.run(start=frame, stop=frame + 1)
    return hba.results.hbonds[:, 2:4].astype(int)


def evaluated_vertex_count(obj: bpy.types.Object) -> int:
    depsgraph = bpy.context.evaluated_depsgraph_get()
    return len(obj.evaluated_get(depsgraph).data.vertices)


def test_hbonds_trajectory(universe):
    mol = mn.Molecule(universe)
    obj = mol.calculations.hbonds().object
    assert mol.frame_attributes["hbonds"] == "CALCULATION"
    assert mol.object.mn.calculations["hbonds"].object == obj
    assert obj.name == f".data_{mol.name}_hbonds"
    assert len(obj.data.vertices) == len(universe.atoms)
    # an independent universe so running the analysis doesn't move the other
    u = mda.Universe(
        data_dir / "md_ppr/box.gro", data_dir / "md_ppr/first_5_frames.xtc"
    )
    for frame in [1, 4, 2]:
        mol.set_frame(frame)
        # the analysis rewinds the trajectory, which must be restored
        assert universe.trajectory.frame == frame
        expected = expected_pairs(
            u,
            frame,
            donors_sel="name N* O*",
            hydrogens_sel="name H*",
            acceptors_sel="name N* O*",
        )
        assert len(expected) > 0
        np.testing.assert_array_equal(edges(obj), expected)


def test_hbonds_structure(fetch):
    # NMR ensemble with hydrogens, bonds and elements
    mol = fetch("2M6Q")
    obj = mol.calculations.hbonds(selection="protein").object
    u = mol.universe.copy()
    for frame in [0, 3]:
        mol.set_frame(frame)
        expected = expected_pairs(
            u,
            frame,
            donors_sel="protein and element N O",
            hydrogens_sel="protein and element H D",
            acceptors_sel="protein and element N O",
        )
        assert len(expected) > 0
        np.testing.assert_array_equal(edges(obj), expected)


def test_hbonds_no_hydrogens(fetch):
    mol = fetch("1l58")
    with pytest.raises(ValueError, match="No hydrogens found"):
        mol.calculations.hbonds()
    assert "hbonds" not in mol.calculations


def test_hbonds_not_updating(fetch):
    mol = fetch("2M6Q")
    obj = mol.calculations.hbonds(updating=False).object
    pairs = edges(obj)
    mol.set_frame(3)
    np.testing.assert_array_equal(edges(obj), pairs)
    # changing a setting recalculates straight away
    mol.calculations.hbonds(d_a_cutoff=3.5)
    assert not np.array_equal(edges(obj), pairs)
    assert mol.calculations._hbonds._kwargs == {"d_a_cutoff": 3.5}
    mol.calculations.hbonds(updating=True)
    mol.set_frame(0)
    u = mol.universe.copy()
    np.testing.assert_array_equal(
        edges(obj),
        expected_pairs(
            u,
            0,
            d_a_cutoff=3.5,
            donors_sel="(all) and (element N O)",
            hydrogens_sel="(all) and (element H D)",
            acceptors_sel="(all) and (element N O)",
        ),
    )


def test_hbonds_remove(fetch):
    mol = fetch("2M6Q")
    obj = mol.calculations.hbonds().object
    n = len(obj.data.edges)
    del mol.calculations["hbonds"]
    assert "hbonds" not in mol.frame_attributes
    assert not mol.calculations._hbonds.is_setup
    # the data object keeps the bonds it last had
    mol.set_frame(3)
    assert len(obj.data.edges) == n


def test_hbonds_pickle(universe):
    mol = mn.Molecule(universe, name="TestPickleHBonds")
    mol.calculations.hbonds("resid 1-50", updating=False, d_a_cutoff=3.2)
    restored = pickle.loads(pickle.dumps(mol))
    # the function isn't pickled, the saved settings are
    assert "hbonds" not in restored.calculations
    restored.set_frame(1)
    hbonds = restored.calculations._hbonds
    assert hbonds.is_setup
    assert hbonds._kwargs == {"d_a_cutoff": 3.2}
    assert restored.props.hbonds.selection == "resid 1-50"
    assert not restored.calculations.updating("hbonds")
    assert hbonds.object == mol.calculations._hbonds.object


def test_style_hbonds_follows_atoms(fetch):
    mol = fetch("2M6Q")
    obj = mol.calculations.hbonds().object
    pairs = edges(obj)
    offset = np.array([1.0, 2.0, 3.0])
    with mol.tree.reset() as (atoms, join):
        (
            atoms
            >> g.SetPosition(offset=tuple(offset))
            >> StyleHBonds(data_object=obj, dashes=1, quality=2)
            >> join
        )
    position = mol.named_attribute("position", evaluate=True)
    # quality 2 gives cylinders with two rings of 8 vertices, one per hydrogen bond
    cylinders = position.reshape(len(pairs), 16, 3)
    atom_position = mol.position + offset
    midpoints = atom_position[pairs].mean(axis=1)
    np.testing.assert_allclose(cylinders.mean(axis=1), midpoints, atol=1e-5)
    # dashes are evenly split along the bond
    with mol.tree.reset() as (atoms, join):
        atoms >> StyleHBonds(data_object=obj, dashes=3, quality=2) >> join
    assert evaluated_vertex_count(mol.object) == len(pairs) * 3 * 16


def test_style_hbonds_selection(fetch):
    mol = fetch("2M6Q")
    obj = mol.calculations.hbonds().object
    pairs = edges(obj)
    is_backbone = mol.named_attribute("is_backbone")
    n_expected = np.all(is_backbone[pairs], axis=1).sum()
    assert 0 < n_expected < len(pairs)
    with mol.tree.reset() as (atoms, join):
        (
            atoms
            >> StyleHBonds(
                selection=g.NamedAttribute.boolean("is_backbone"),
                data_object=obj,
                dashes=1,
                quality=2,
            )
            >> join
        )
    assert evaluated_vertex_count(mol.object) == n_expected * 16


def test_hbonds_init_operator(fetch):
    mol = fetch("2M6Q")
    hbonds = mol.calculations._hbonds
    assert not hbonds.is_setup
    bpy.ops.mn.hbonds_init(uuid=mol.uuid, selection="protein")
    assert hbonds.is_setup
    assert mol.props.hbonds.selection == "protein"
    n = len(hbonds.object.data.edges)
    # a second run updates the selection
    bpy.ops.mn.hbonds_init(uuid=mol.uuid, selection="resid 1-20")
    assert mol.props.hbonds.selection == "resid 1-20"
    assert len(hbonds.object.data.edges) < n
