from typing import cast
import bpy
import MDAnalysis as mda
import numpy as np
import pytest
from MDAnalysis.lib.distances import apply_PBC
from MDAnalysis.lib.mdamath import triclinic_vectors
from MDAnalysis.tests.datafiles import DCD, GRO, PSF, XTC
from nodebpy.nodes.geometry import SetPosition, StoreNamedAttribute
import molecularnodes as mn
from molecularnodes.nodes.geometry import FragmentID, PeriodicOffset, PeriodicUnwrap
from .constants import data_dir

BOXES = {
    "cubic": [85.5, 85.5, 85.5, 90.0, 90.0, 90.0],
    "orthorhombic": [53.8, 44.0, 57.2, 90.0, 90.0, 90.0],
    "hexagonal": [98.0, 98.0, 130.7, 90.0, 90.0, 120.0],
    "dodecahedron": [82.68, 82.68, 82.68, 60.0, 60.0, 90.0],
}
BOX_NAMES = ["a", "b", "c", "alpha", "beta", "gamma"]
# Å, and the tolerance for comparing positions read back from the mesh
ATOL = 1e-2
# fractional offset that splits molecules across the box before unwrapping
SHIFT = np.array([0.37, 0.21, 0.5])


def _box_inputs(box) -> dict:
    return {"update": False, **dict(zip(BOX_NAMES, (float(x) for x in box)))}


def _box_values(node) -> list[float]:
    """The box lengths and angles set on a node's inputs"""
    return [
        cast(bpy.types.NodeSocketFloat, node.node.inputs[name]).default_value
        for name in BOX_NAMES
    ]


def _evaluated_positions(mol: mn.Molecule) -> np.ndarray:
    """Positions after the node tree, in Å"""
    return mol.named_attribute("position", evaluate=True) / mol.world_scale


def _assert_inside_box(positions: np.ndarray, box) -> None:
    fractional = positions @ np.linalg.inv(triclinic_vectors(box))
    assert fractional.min() > -1e-4 and fractional.max() < 1 + 1e-4


def _assert_same_image(result: np.ndarray, expected: np.ndarray, box) -> None:
    """Equal, allowing for atoms on a box face wrapped to the opposite one"""
    fractional = (result - expected) @ np.linalg.inv(triclinic_vectors(box))
    fractional -= np.round(fractional)
    assert np.abs(fractional @ triclinic_vectors(box)).max() < ATOL


def _assert_same_up_to_lattice(result, expected, groups, box) -> None:
    """Each group is the same shape as expected, moved by a whole box vector"""
    _, first, inverse = np.unique(groups, return_index=True, return_inverse=True)
    diff = result - expected
    shape_error = diff - diff[first][inverse]
    assert np.abs(shape_error).max() < ATOL
    lattice = diff[first] @ np.linalg.inv(triclinic_vectors(box))
    assert np.allclose(lattice, np.round(lattice), atol=1e-3)


@pytest.fixture
def mol_triclinic():
    # rhombic dodecahedron with bonds, whole molecules and 11779 fragments
    u = mda.Universe(data_dir / "md_ppr/md.tpr", data_dir / "md_ppr/md.gro")
    return mn.Molecule(u)


@pytest.mark.parametrize("box_name", BOXES)
@pytest.mark.parametrize(
    "offset_type, offset",
    [
        ("Factor", (0.0, 0.0, 0.0)),
        ("Factor", (0.37, -0.21, 1.5)),
        ("Distance", (1.3, -2.2, 0.4)),
    ],
)
def test_periodic_offset(mol_triclinic, box_name, offset_type, offset):
    mol = mol_triclinic
    box = np.array(BOXES[box_name])
    positions = mol.named_attribute("position") / mol.world_scale
    with mol.tree.reset() as (atoms, join):
        node = PeriodicOffset(
            offset_type=offset_type,
            offset_factor=offset,
            offset_distance=offset,
            **_box_inputs(box),
        )
        atoms >> SetPosition(position=node) >> join
    if offset_type == "Factor":
        shift = np.array(offset) @ triclinic_vectors(box)
    else:
        # world units, where 1 unit is 1 nm
        shift = np.array(offset) / mol.world_scale
    expected = apply_PBC((positions + shift).astype(np.float32), box)
    result = _evaluated_positions(mol)
    _assert_inside_box(result, box)
    _assert_same_image(result, expected, box)


@pytest.mark.parametrize("node", [PeriodicOffset, PeriodicUnwrap])
def test_periodic_no_box_unchanged(mol_triclinic, node):
    mol = mol_triclinic
    positions = mol.named_attribute("position")
    with mol.tree.reset() as (atoms, join):
        atoms >> SetPosition(position=node(update=False)) >> join
    np.testing.assert_allclose(
        mol.named_attribute("position", evaluate=True), positions, atol=1e-6
    )


def test_fragment_id(mol_triclinic):
    mol = mol_triclinic
    expected = mol.universe.atoms.fragindices
    np.testing.assert_array_equal(mol.named_attribute("fragment_id"), expected)
    # read by the Fragment ID node
    with mol.tree.reset() as (atoms, join):
        store = StoreNamedAttribute.point.integer(
            name="node_fragment_id", value=FragmentID()
        )
        atoms >> store >> join
    np.testing.assert_array_equal(
        mol.named_attribute("node_fragment_id", evaluate=True), expected
    )
    # a topology without bonds has no fragments
    assert "fragment_id" not in mn.Molecule(mda.Universe(GRO)).list_attributes()


def _break_and_unwrap(mol, box, place_groups, group_id=lambda: None):
    """Split molecules across the box with Periodic Offset, then unwrap them"""
    with mol.tree.reset() as (atoms, join):
        offset = PeriodicOffset(offset_factor=tuple(SHIFT), **_box_inputs(box))
        unwrap = PeriodicUnwrap(
            group_id=group_id(), place_groups=place_groups, **_box_inputs(box)
        )
        (atoms >> SetPosition(position=offset) >> SetPosition(position=unwrap) >> join)
    positions = mol.named_attribute("position") / mol.world_scale
    broken = apply_PBC(
        (positions + SHIFT @ triclinic_vectors(box)).astype(np.float32),
        box,
    )
    return broken, _evaluated_positions(mol)


@pytest.mark.parametrize("place_groups", [False, True])
def test_periodic_unwrap(mol_triclinic, place_groups):
    mol = mol_triclinic
    u = mol.universe
    box = u.dimensions
    positions = mol.named_attribute("position") / mol.world_scale
    broken, result = _break_and_unwrap(mol, box, place_groups, FragmentID)
    # the offset splits molecules, which the unwrap puts back together
    assert u.bonds.values().max() < 2.0
    u.atoms.positions = broken
    assert u.bonds.values().max() > 10.0
    u.atoms.positions = result
    assert u.bonds.values().max() < 2.0
    fragments = u.atoms.fragindices
    shifted = positions + SHIFT @ triclinic_vectors(box)
    _assert_same_up_to_lattice(result, shifted, fragments, box)
    if place_groups:
        u.atoms.positions = result
        centres = u.atoms.center_of_geometry(compound="fragments")
        fractional = centres @ np.linalg.inv(triclinic_vectors(box))
        assert fractional.min() > -0.1 and fractional.max() < 1.1


@pytest.mark.parametrize("place_groups", [False, True])
def test_periodic_unwrap_default_group_per_atom(mol_triclinic, place_groups):
    # unlinked, Group ID is the ID or index, so each atom is its own group
    mol = mol_triclinic
    box = mol.universe.dimensions
    broken, result = _break_and_unwrap(mol, box, place_groups)
    _assert_same_image(result, broken, box)
    if place_groups:
        _assert_inside_box(result, box)
    else:
        assert np.abs(result - broken).max() < ATOL


def test_periodic_box_updates():
    mol = mn.Molecule.load(GRO, XTC)
    with mol.tree.reset() as (atoms, join):
        offset = PeriodicOffset()
        unwrap = PeriodicUnwrap()
        atoms >> SetPosition(position=offset) >> SetPosition(position=unwrap) >> join
    for frame in [1, 5]:
        mol.set_frame(frame)
        dimensions = mol.universe.trajectory.ts.dimensions
        for node in [offset, unwrap]:
            np.testing.assert_allclose(_box_values(node), dimensions, rtol=1e-6)


def test_periodic_box_updates_without_dimensions():
    mol = mn.Molecule.load(PSF, DCD)
    with mol.tree.reset() as (atoms, join):
        node = PeriodicOffset()
        atoms >> SetPosition(position=node) >> join
    mol.set_frame(1)
    assert _box_values(node) == [0] * 6


def _dataset_universe(name: str) -> mda.Universe:
    """A universe from MDAnalysisData, skipping when it isn't downloaded"""
    from MDAnalysisData import datasets

    try:
        data = getattr(datasets, f"fetch_{name}")(download_if_missing=False)
    except OSError:
        pytest.skip(f"MDAnalysisData '{name}' is not downloaded")
    return mda.Universe(data.topology, data.trajectory)


@pytest.mark.parametrize(
    "name",
    [
        "adk_equilibrium",
        "PEG_1chain",
        "membrane_peptide",
        "yiip_equilibrium_short",
        "ifabp_water",
    ],
)
@pytest.mark.parametrize("place_groups", [False, True])
def test_periodic_unwrap_datasets(name, place_groups):
    u = _dataset_universe(name)
    if name == "yiip_equilibrium_short":
        u.guess_TopologyAttrs(to_guess=["bonds"], vdwradii={"Z": 1.39, "D": 1.2})
    mol = mn.Molecule(u)
    box = u.dimensions
    if box is None:
        # no box, as in the fitted ifabp_water trajectory: nothing changes
        positions = mol.named_attribute("position")
        with mol.tree.reset() as (atoms, join):
            node = PeriodicUnwrap(group_id=FragmentID())
            atoms >> SetPosition(position=node) >> join
        np.testing.assert_allclose(
            mol.named_attribute("position", evaluate=True), positions, atol=1e-6
        )
        return
    broken, result = _break_and_unwrap(mol, box, place_groups, FragmentID)
    atoms = u.atoms
    assert atoms is not None
    atoms.positions = broken
    expected = atoms.unwrap(compound="fragments")
    _assert_same_up_to_lattice(result, expected, atoms.fragindices, box)
