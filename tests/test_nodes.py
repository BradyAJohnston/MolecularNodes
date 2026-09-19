import random
from typing import Any
import bpy
import MDAnalysis as mda
import numpy as np
import pytest
from databpy.nodes import get_input, get_output
from MDAnalysis.tests.datafiles import DCD, GRO, PSF, XTC
from nodebpy.nodes.geometry import (
    Compare,
    GetBundleItem,
    GetGeometryBundle,
    Group,
    Index,
    MeshLine,
    Points,
    RealizeInstances,
    SetPosition,
    StoreNamedAttribute,
)
import molecularnodes as mn
from molecularnodes.nodes._utils import (
    custom_boolean_iswitch,
    custom_color_iswitch,
    get_final_style_nodes,
)
from molecularnodes.nodes.geometry import (
    BreakBonds,
    BuildElasticNetwork,
    Charge,
    ColorToOKLab,
    FindBonds,
    NucleicChi,
    NucleicDihedral,
    OKLabToColor,
    PeptideChi,
    PeptideDihedral,
    PeriodicArray,
    SegmentID,
    SetColor,
    SimulateElasticNetwork,
    StyleCartoon,
)
from .constants import codes, data_dir
from .utils import GeometrySet, NumpySnapshotExtension

random.seed(6)


def _input_sockets(tree: bpy.types.NodeTree) -> dict:
    return {
        item.name: item
        for item in tree.interface.items_tree
        if item.item_type == "SOCKET" and item.in_out == "INPUT"
    }


def test_get_nodes():
    mol = mn.Molecule.fetch("4ozs", cache=data_dir).add_style("spheres")

    last = get_output(mol.node_group).inputs[0].links[0].from_node
    assert last.name == "Join Geometry"
    style = get_final_style_nodes(mol.node_group)[0]
    assert style.name == "Style Spheres"
    assert style.node_tree.name == "Style Spheres"

    mol2 = mn.Molecule.fetch("1cd3", cache=data_dir).add_style("cartoon")

    assert get_final_style_nodes(mol2.node_group)[0].node_tree.name == "Style Cartoon"


def test_selection():
    chain_ids = [let for let in "ABCDEFG123456"]
    tree = custom_boolean_iswitch("test_node", chain_ids, prefix="Chain ")

    input_sockets = _input_sockets(tree.tree)
    for letter, socket in zip(chain_ids, input_sockets.values()):
        assert f"Chain {letter}" == socket.name
        assert socket.default_value is False


@pytest.mark.parametrize("code", codes)
@pytest.mark.parametrize("attribute", ["chain_id", "entity_id"])
def test_selection_working(snapshot_custom: NumpySnapshotExtension, attribute, code):
    mol = mn.Molecule.fetch(code, cache=data_dir)
    with mol.tree.reset() as (atoms, join):
        sel = Group()
        sel.node.node_tree = custom_boolean_iswitch(
            mol.name, getattr(mol.props, f"{attribute}s"), attribute
        ).tree
        atoms >> StyleCartoon(selection=sel) >> RealizeInstances() >> join

    _n = len(sel.i)

    for inp in sel.i:
        inp.default_value = True  # type: ignore
        pos = mol.named_attribute("position", evaluate=True)
        assert snapshot_custom == pos.shape
        assert snapshot_custom == pos
        inp.default_value = False


@pytest.mark.parametrize("code", codes)
@pytest.mark.parametrize("attribute", ["chain_id", "entity_id"])
def test_color_custom(snapshot_custom: NumpySnapshotExtension, code, attribute):
    mol = mn.Molecule.fetch(code, cache=data_dir)

    group_col = custom_color_iswitch(
        name=f"Color Entity {mol.name}",
        items=getattr(mol.props, f"{attribute}s"),
        attribute_name=attribute,
    )
    with mol.tree.reset() as (atoms, join):
        n_color = Group()
        n_color.node.node_tree = group_col.tree

        atoms >> SetColor(color=n_color) >> StyleCartoon() >> join

    for i, input in enumerate(n_color.i):
        setattr(input, "default_value", mn.color.random_rgb(i))

    assert snapshot_custom == mol.named_attribute("Color")


def test_iswitch_creation():
    items = [str(x) for x in range(10)]
    tree_boolean = custom_boolean_iswitch("newboolean", items).tree
    # ensure there isn't an item called 'Color' in the created interface
    assert not tree_boolean.interface.items_tree.get("Color")
    assert tree_boolean.interface.items_tree["Selection"].in_out == "OUTPUT"
    assert tree_boolean.interface.items_tree["Inverted"].in_out == "OUTPUT"
    for i in items:
        assert tree_boolean.interface.items_tree[str(i)].in_out == "INPUT"

    tree_rgba = custom_color_iswitch("newcolor", items).tree
    # ensure there isn't an item called 'selection'
    assert not tree_rgba.interface.items_tree.get("Selection")
    assert tree_rgba.interface.items_tree["Color"].in_out == "OUTPUT"
    for i in items:
        assert tree_rgba.interface.items_tree[str(i)].in_out == "INPUT"


def test_op_custom_color():
    mol = mn.Molecule.load(data_dir / "1cd3.cif")
    mol.object.select_set(True)
    group = custom_color_iswitch(
        name=f"Color Chain {mol.name}", items=mol.props.chain_ids
    ).tree

    assert group
    assert group.interface.items_tree["G"].name == "G"
    assert group.interface.items_tree[-1].name == "G"
    assert group.interface.items_tree[0].name == "Color"


def test_color_lookup_supplied():
    col = mn.color.random_rgb(6)
    name = "test"
    node = custom_color_iswitch(
        name=name,
        items={str(x): col for x in range(10, 20)},
        offset=10,
    ).tree
    assert node.name == name
    for item in _input_sockets(node).values():
        assert np.allclose(np.array(item.default_value), col)

    node = custom_color_iswitch(name="test2", items=range(10, 20), offset=10).tree
    for item in _input_sockets(node).values():
        assert not np.allclose(np.array(item.default_value), col)


@pytest.mark.parametrize(
    "node", [PeptideDihedral, NucleicDihedral, PeptideChi, NucleicChi]
)
@pytest.mark.parametrize("code", ["8H1B", "1BNA"])
def test_dihedral_rotations(snapshot_custom: NumpySnapshotExtension, code, node):
    mol = mn.Molecule.fetch(code, cache=data_dir)
    with mol.tree.reset() as (atoms, join):
        pos = node()
        (atoms >> SetPosition(position=pos) >> join)

    for input in pos.i:
        if input.name in ["Position", "Selection"]:
            continue
        input.default_value = 1.0
    assert snapshot_custom == mol.named_attribute("position", evaluate=True)[:100]


def test_topo_bonds():
    mol = mn.Molecule.fetch("1BNA", cache=data_dir)
    with mol.tree.reset() as (atoms, join):
        atoms >> BreakBonds(cutoff=0.0) >> join

    # compare the number of edges before and after deleting them with
    gs = GeometrySet(mol.object)
    assert gs.mesh
    assert len(gs.mesh.edges) == 0

    # add the node to find the bonds, and ensure the number of bonds pre and post the nodes
    # are the same (other attributes will be different, but for now this is good)
    with mol.tree.reset() as (atoms, join):
        atoms >> BreakBonds(cutoff=0.0) >> FindBonds() >> join

    gs_new = GeometrySet(mol.object)
    assert len(mol.object.data.edges) == len(gs_new.mesh.edges)


def test_is_modifier():
    bpy.ops.wm.open_mainfile(filepath=str(mn.assets.MN_DATA_FILE))
    for tree in bpy.data.node_groups:
        if tree.name.startswith("Style") and "Preset" not in tree.name:
            assert tree.is_modifier
    mol = mn.Molecule.fetch("4ozs").add_style("spheres")
    assert mol.modifier_node_tree.is_modifier


def test_node_setup():
    mn.Molecule.fetch("4ozs").add_style("spheres")
    tree = bpy.data.node_groups["MN_4ozs"]
    assert tree.interface.items_tree["Atoms"].name == "Atoms"
    assert list(get_input(tree).outputs.keys()) == ["Atoms", ""]
    assert list(get_output(tree).inputs.keys()) == ["Geometry", ""]


def test_reuse_node_group():
    mol = mn.Molecule.fetch("4ozs").add_style("spheres")
    tree = bpy.data.node_groups["MN_4ozs"]
    n_nodes = len(tree.nodes)
    bpy.data.objects.remove(mol.object)
    del mol
    assert n_nodes == len(tree.nodes)
    mn.Molecule.fetch("4ozs")
    assert n_nodes == len(tree.nodes)


def _get_node_defaults(node) -> list[Any]:
    defaults = []
    for input in node.inputs:
        if not hasattr(input, "default_value"):
            continue
        default = input.default_value
        if isinstance(default, float):
            default = round(default, 3)
        defaults.append(default)

    return defaults


def test_periodic_array(snapshot, tmp_path):
    traj = mn.Molecule.load(GRO, XTC)

    with traj.tree.reset() as (atoms, join):
        node = PeriodicArray()
        atoms >> node >> join

    traj.set_frame(1)
    defaults_0 = _get_node_defaults(node.node)
    traj.set_frame(10)
    defaults_10 = _get_node_defaults(node.node)

    dim_idx = slice(1, 7)
    assert not all([x == y for x, y in zip(defaults_0[dim_idx], defaults_10[dim_idx])])
    # the unit cell dimensions ar currently inputs 1..7 for the node as it is setup so
    # we just subset those and check it matches the universe
    assert np.allclose(defaults_10[dim_idx], traj.universe.trajectory.ts.dimensions)

    # for some reason we need to trigger a proper re-evaluation of the GN node tree
    # by saving to a temp file #TODO: look into and try to fix this
    bpy.ops.wm.save_as_mainfile(filepath=str(tmp_path / "example.blend"))
    assert snapshot == GeometrySet(traj.object).summary()


# this topology doesn't have any dimension information so it should just
# update the positions and _attempt_ to update the periodic box but fail not do so quietly
# and everything remains 0
def test_periodic_array_no_dimensions():
    traj = mn.Molecule.load(PSF, DCD)
    with traj.tree.reset() as (atoms, join):
        node = PeriodicArray()
        atoms >> node >> join

    traj.set_frame(1)
    defaults_0 = _get_node_defaults(node.node)
    traj.set_frame(frame=10)
    defaults_10 = _get_node_defaults(node.node)

    assert defaults_0 == defaults_10
    assert defaults_0[1:7] == [0] * 6


def _store_charge_node(mol):
    with mol.tree.reset() as (atoms, join):
        (
            atoms
            >> StoreNamedAttribute.point.float(name="charge_node", value=Charge())
            >> join
        )
    return mol.named_attribute("charge_node", evaluate=True)


def test_charge_node():
    # 4ozs provides no charges, so there is no `charge` attribute and the node reads 0.0
    mol = mn.Molecule.fetch("4ozs", cache=data_dir)
    assert "charge" not in mol.list_attributes()
    assert np.all(_store_charge_node(mol) == 0)

    # 8U8W carries formal charges on its ions, which the node reads back verbatim
    mol = mn.Molecule.fetch("8U8W", cache=data_dir)
    expected = mol.named_attribute("charge")
    assert np.any(expected != 0)
    assert np.allclose(_store_charge_node(mol), expected)


@pytest.mark.parametrize(
    "topology, trajectory, n_segments",
    [
        ("md_ppr/md.tpr", "md_ppr/md.gro", 3),
        ("md_ppr/box.gro", "md_ppr/first_5_frames.xtc", 1),
    ],
)
def test_segment_id(topology, trajectory, n_segments):
    universe = mda.Universe(data_dir / topology, data_dir / trajectory)
    traj = mn.Molecule(universe)
    expected = traj.named_attribute("segid")
    assert len(np.unique(expected)) == n_segments
    assert np.array_equal(np.unique(expected), np.arange(n_segments))

    with traj.tree.reset() as (atoms, join):
        (
            atoms
            >> StoreNamedAttribute.point.integer(name="node_segid", value=SegmentID())
            >> join
        )

    assert np.array_equal(traj.named_attribute("node_segid", evaluate=True), expected)


def test_build_elastic_network():
    mol = mn.Molecule.fetch("4ozs")

    with mol.tree.reset() as (atoms, join):
        BuildElasticNetwork(atoms) >> join

    gs = GeometrySet(mol.object)
    assert gs.mesh
    assert len(gs.mesh.edges) == 1049
    assert len(gs.mesh.vertices) == sum(mol["is_alpha_carbon"])


def test_evaluate_on_atoms_bundle():
    """
    The `Evaluate on Atoms` node was previously still storing the 'MN/Atoms' bundle,
    even when the output was set to geometry. It was just filling it with emptry geometry.
    This was leading to the style nodes failing to work properly if we were finding bonds,
    as the resulting bonded geoemtry was output as the `Geometry` but the 'MN/Atoms'
    bundle was populated with the original pre-bonded geometry, meaning the style used
    this old / outdated information instead of the results of the bond calculation.

    For the test we want to check that the bundle isn't being created fromt he `FindBonds()`
    and we check if it exists and create a point if so which we can test for with pytest.
    """

    mol = mn.Molecule.fetch("4ozs")

    with mol.tree.reset() as (atoms, join):
        bundle = (atoms >> FindBonds() >> GetGeometryBundle()).o.bundle

        # importantly we have to use the typed "GetBundleItem" becuase
        # if the item we are requesting doesn't have the same type the "exists"
        # returns false, even if it exists but is of a different type
        Points(GetBundleItem.bundle(bundle, "MN").o.exists) >> join

    gs = GeometrySet(mol.object)
    assert gs.pointcloud is None or len(gs.pointcloud.points) == 0


# Reference values from Ottosson, "A perceptual color space for image processing"
# (https://bottosson.github.io/posts/oklab/), linear sRGB in, OKLab out.
OKLAB_REFERENCE = {
    (1.0, 0.0, 0.0): (0.62796, 0.22486, 0.12585),
    (0.0, 1.0, 0.0): (0.86644, -0.23389, 0.17950),
    (0.0, 0.0, 1.0): (0.45201, -0.03246, -0.31153),
    (1.0, 1.0, 1.0): (1.0, 0.0, 0.0),
}


@pytest.mark.parametrize("rgb", list(OKLAB_REFERENCE))
def test_color_to_oklab(rgb):
    """Color to OKLab matches Ottosson's reference values and round-trips."""
    mol = mn.Molecule.fetch("4ozs", cache=data_dir)
    with mol.tree.reset() as (atoms, join):
        oklab = ColorToOKLab(color=(*rgb, 1.0))
        (
            atoms
            >> StoreNamedAttribute.point.vector(name="oklab", value=oklab)
            >> StoreNamedAttribute.point.color(
                name="rgb", value=OKLabToColor(oklab=oklab)
            )
            >> join
        )
    lab = mol.named_attribute("oklab", evaluate=True)[0]
    back = mol.named_attribute("rgb", evaluate=True)[0, :3]
    assert np.allclose(lab, OKLAB_REFERENCE[rgb], atol=1e-3), lab
    assert np.allclose(back, rgb, atol=1e-4), back


def _simulate_two_point_network(mol, frames: int):
    """Two points 1.0 apart with masses 1 and 3, joined by one edge whose rest
    length is set to 0.5, simulated with no external forces."""
    with mol.tree.reset() as (atoms, join):
        mass = Compare.integer.equal(Index().o.index, 0).o.result.switch.float(3.0, 1.0)
        (
            MeshLine(count=2, start_location=(0.0, 0.0, 0.0), offset=(1.0, 0.0, 0.0))
            >> StoreNamedAttribute.point.float(name="mass", value=mass)
            >> SimulateElasticNetwork(
                substeps=1,
                force=(0.0, 0.0, 0.0),
                drag=0.0,
                edge_length_source="Custom",
                edge_length=0.5,
            )
            >> join
        )
    scene = bpy.context.scene
    start = scene.frame_current
    try:
        for f in range(start, start + frames + 1):
            scene.frame_set(f)
        return (
            mol.named_attribute("position", evaluate=True),
            mol.named_attribute("inverse_mass", evaluate=True),
        )
    finally:
        scene.frame_set(start)


def test_simulate_elastic_network_inverse_mass():
    """The stored inverse mass is 1 / mass."""
    mol = mn.Molecule.fetch("4ozs", cache=data_dir)
    _, inverse_mass = _simulate_two_point_network(mol, frames=0)
    assert np.allclose(inverse_mass, [1.0, 1.0 / 3.0], atol=1e-6), inverse_mass


def test_simulate_elastic_network_two_points():
    """One XPBD step of a single distance constraint: each point moves in
    proportion to its inverse mass, the mass-weighted centre stays put and the
    edge reaches its rest length."""
    mol = mn.Molecule.fetch("4ozs", cache=data_dir)
    positions, _ = _simulate_two_point_network(mol, frames=1)
    # C = 1.0 - 0.5, point 0 (w=1) moves C * 1 / (1 + 1/3), point 1 (w=1/3) moves
    # C * (1/3) / (1 + 1/3), both toward each other along x
    expected = np.array([[0.375, 0.0, 0.0], [0.875, 0.0, 0.0]])
    assert np.allclose(positions, expected, atol=1e-4), positions
    masses = np.array([1.0, 3.0])
    centre = (positions * masses[:, None]).sum(axis=0) / masses.sum()
    assert np.allclose(centre, [0.75, 0.0, 0.0], atol=1e-4), centre
    assert np.isclose(np.linalg.norm(positions[1] - positions[0]), 0.5, atol=1e-4)
