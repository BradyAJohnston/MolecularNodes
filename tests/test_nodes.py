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
    NamedAttribute,
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
    AnimateAction,
    BreakBonds,
    BuildElasticNetwork,
    Charge,
    ColorToOKLab,
    EaseValue,
    FadeGeometry,
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
    StaggerValue,
    StyleCartoon,
    SymmetryCyclic,
    SymmetryDihedral,
    SymmetryHelical,
    UResID,
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


# Generated symmetry that reproduces each entry's deposited biological assembly:
# the parameters were fitted from the entry's pdbx_struct_oper_list operators.
# Shared with the golden renders in test_render_images.py.
SYMMETRY_EXAMPLES = {
    # streptavidin, D2 tetramer built from one chain
    "1STP": lambda: SymmetryDihedral(
        order=2, axis=(-0.707107, 0.707107, 0.0), offset=-3 * np.pi / 2
    ),
    # DatZ phosphohydrolase, D3 hexamer from one chain
    "6ZPA": lambda: SymmetryDihedral(order=3, offset=-np.pi / 3),
    # HIV-1 capsid hexamer, C6 from one chain
    "8QUK": lambda: SymmetryCyclic(order=6),
    # E. coli Hfq, C6 ring whose six-fold axis is off the origin
    "4PNO": lambda: SymmetryCyclic(order=6, centre=(6.111, 0.0, 0.0)),
    # tobacco mosaic virus, 49-subunit helix: 1.408 A rise, 22.03 degree twist
    "4UDV": lambda: SymmetryHelical(
        count=49, rise=1.408, twist=0.384496, axis=(0.0, 0.0, -1.0)
    ),
}


# the helix is deposited as k = -24..24 about the model, which the node builds
# as k = 0..48 instead, so it has no deposited assembly to compare against
@pytest.mark.parametrize("code", ["1STP", "6ZPA", "8QUK", "4PNO"])
def test_symmetry_matches_deposited_assembly(code):
    from scipy.spatial import cKDTree

    mol = mn.Molecule.fetch(code)
    positions = mol.named_attribute("position")
    deposited = np.vstack(
        [
            positions @ np.array(op["matrix"])[:3, :3].T
            + np.array(op["matrix"])[:3, 3] * mol.world_scale
            for op in mol.assemblies()["1"]
        ]
    )

    with mol.tree.reset() as (atoms, join):
        atoms >> SYMMETRY_EXAMPLES[code]() >> RealizeInstances() >> join
    generated = mol.named_attribute("position", evaluate=True)

    assert len(generated) == len(deposited)
    distance, _ = cKDTree(deposited).query(generated)
    assert distance.max() < 1e-3


def _rotation_about(axis, angle: float) -> np.ndarray:
    "A 3x3 rotation of `angle` radians about `axis` (Rodrigues)."
    a = np.asarray(axis, dtype=float)
    a = a / np.linalg.norm(a)
    cross = np.array([[0.0, -a[2], a[1]], [a[2], 0.0, -a[0]], [-a[1], a[0], 0.0]])
    return np.eye(3) + np.sin(angle) * cross + (1.0 - np.cos(angle)) * (cross @ cross)


def _expected_copies(positions, operators, centre) -> np.ndarray:
    "Apply each (rotation, translation) about `centre`, concatenated in order."
    centre = np.asarray(centre, dtype=float)
    return np.vstack(
        [
            (positions - centre) @ rotation.T + centre + translation
            for rotation, translation in operators
        ]
    )


def _realized_symmetry(make_node) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    "Original and realized positions, realized sym_id and world scale for a symmetry node."
    mol = mn.Molecule.fetch("4ozs", cache=data_dir)
    positions = mol.named_attribute("position")
    with mol.tree.reset() as (atoms, join):
        atoms >> make_node() >> RealizeInstances() >> join
    return (
        positions,
        mol.named_attribute("position", evaluate=True),
        mol.named_attribute("sym_id", evaluate=True),
        mol.world_scale,
    )


def test_symmetry_cyclic():
    axis, centre = (0.3, -0.2, 1.0), (0.5, 0.1, -0.2)
    positions, realized, sym_id, _ = _realized_symmetry(
        lambda: SymmetryCyclic(order=5, axis=axis, centre=centre)
    )
    operators = [
        (_rotation_about(axis, 2 * np.pi * k / 5), np.zeros(3)) for k in range(5)
    ]
    assert realized.shape == (len(positions) * 5, 3)
    assert np.allclose(
        realized, _expected_copies(positions, operators, centre), atol=1e-4
    )
    assert np.array_equal(sym_id, np.repeat(np.arange(5), len(positions)))


def test_symmetry_cyclic_factor():
    axis, centre = (0.0, 0.0, 1.0), (0.5, 0.1, -0.2)
    positions, realized, _, _ = _realized_symmetry(
        lambda: SymmetryCyclic(order=4, axis=axis, centre=centre, factor=0.0)
    )
    assert np.allclose(realized, np.tile(positions, (4, 1)), atol=1e-4)

    positions, realized, _, _ = _realized_symmetry(
        lambda: SymmetryCyclic(order=4, axis=axis, centre=centre, factor=0.5)
    )
    operators = [
        (_rotation_about(axis, 0.5 * 2 * np.pi * k / 4), np.zeros(3)) for k in range(4)
    ]
    assert np.allclose(
        realized, _expected_copies(positions, operators, centre), atol=1e-4
    )


def test_symmetry_dihedral():
    axis, centre = (0.0, 0.0, 1.0), (0.2, -0.4, 0.1)
    positions, realized, sym_id, _ = _realized_symmetry(
        lambda: SymmetryDihedral(order=3, axis=axis, centre=centre)
    )
    # the two-fold for a Z axis lies along Y, matching ProteinBlender's choice
    flip = _rotation_about((0.0, 1.0, 0.0), np.pi)
    ring = [_rotation_about(axis, 2 * np.pi * k / 3) for k in range(3)]
    operators = [(r, np.zeros(3)) for r in ring] + [
        (r @ flip, np.zeros(3)) for r in ring
    ]
    assert realized.shape == (len(positions) * 6, 3)
    assert np.allclose(
        realized, _expected_copies(positions, operators, centre), atol=1e-4
    )
    assert np.array_equal(sym_id, np.repeat(np.arange(6), len(positions)))


def test_symmetry_dihedral_offset():
    # the offset rotates the flipped ring about the axis, moving the two-fold by half of it
    axis, centre, offset = (0.0, 0.0, 1.0), (0.0, 0.0, 0.0), 0.7
    positions, realized, _, _ = _realized_symmetry(
        lambda: SymmetryDihedral(order=2, axis=axis, centre=centre, offset=offset)
    )
    flip = _rotation_about((0.0, 1.0, 0.0), np.pi)
    ring = [_rotation_about(axis, np.pi * k) for k in range(2)]
    operators = [(r, np.zeros(3)) for r in ring] + [
        (_rotation_about(axis, np.pi * k + offset) @ flip, np.zeros(3))
        for k in range(2)
    ]
    assert np.allclose(
        realized, _expected_copies(positions, operators, centre), atol=1e-4
    )


def test_symmetry_helical():
    axis, centre = (0.0, 1.0, 1.0), (0.1, 0.2, 0.3)
    rise, twist = 27.5, np.radians(-166.7)
    positions, realized, sym_id, world_scale = _realized_symmetry(
        lambda: SymmetryHelical(
            count=6, rise=rise, twist=twist, axis=axis, centre=centre
        )
    )
    direction = np.asarray(axis) / np.linalg.norm(axis)
    operators = [
        (_rotation_about(axis, twist * k), direction * rise * world_scale * k)
        for k in range(6)
    ]
    assert realized.shape == (len(positions) * 6, 3)
    assert np.allclose(
        realized, _expected_copies(positions, operators, centre), atol=1e-4
    )
    assert np.array_equal(sym_id, np.repeat(np.arange(6), len(positions)))


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


# Robert Penner's easing functions as listed on https://easings.net, ported
# separately for In, Out and In Out so the node's shared-curve construction
# (Out = 1 - f(1 - t), In Out from f(2t) and f(2 - 2t)) is checked against the
# reference rather than against itself.
def _penner(curve: str, ease: str, x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    pi = np.pi
    c1 = 1.70158
    c2 = c1 * 1.525
    c3 = c1 + 1
    c4 = 2 * pi / 3
    c5 = 2 * pi / 4.5

    def bounce_out(u):
        n1, d1 = 7.5625, 2.75
        return np.where(
            u < 1 / d1,
            n1 * u * u,
            np.where(
                u < 2 / d1,
                n1 * (u - 1.5 / d1) ** 2 + 0.75,
                np.where(
                    u < 2.5 / d1,
                    n1 * (u - 2.25 / d1) ** 2 + 0.9375,
                    n1 * (u - 2.625 / d1) ** 2 + 0.984375,
                ),
            ),
        )

    def poly(n):
        return {
            "In": x**n,
            "Out": 1 - (1 - x) ** n,
            "In Out": np.where(x < 0.5, 2 ** (n - 1) * x**n, 1 - (-2 * x + 2) ** n / 2),
        }

    def pinned(inner, x0=0.0, x1=1.0):
        return np.where(x == 0, x0, np.where(x == 1, x1, inner))

    with np.errstate(invalid="ignore", divide="ignore"):
        table = {
            "Linear": {"In": x, "Out": x, "In Out": x},
            "Sinusoidal": {
                "In": 1 - np.cos(x * pi / 2),
                "Out": np.sin(x * pi / 2),
                "In Out": -(np.cos(pi * x) - 1) / 2,
            },
            "Quadratic": poly(2),
            "Cubic": poly(3),
            "Quartic": poly(4),
            "Quintic": poly(5),
            "Exponential": {
                "In": np.where(x == 0, 0.0, 2 ** (10 * x - 10)),
                "Out": np.where(x == 1, 1.0, 1 - 2 ** (-10 * x)),
                "In Out": pinned(
                    np.where(
                        x < 0.5,
                        2 ** (20 * x - 10) / 2,
                        (2 - 2 ** (-20 * x + 10)) / 2,
                    )
                ),
            },
            "Circular": {
                "In": 1 - np.sqrt(1 - x**2),
                "Out": np.sqrt(1 - (x - 1) ** 2),
                "In Out": np.where(
                    x < 0.5,
                    (1 - np.sqrt(1 - (2 * x) ** 2)) / 2,
                    (np.sqrt(1 - (-2 * x + 2) ** 2) + 1) / 2,
                ),
            },
            "Back": {
                "In": c3 * x**3 - c1 * x**2,
                "Out": 1 + c3 * (x - 1) ** 3 + c1 * (x - 1) ** 2,
                "In Out": np.where(
                    x < 0.5,
                    ((2 * x) ** 2 * ((c2 + 1) * 2 * x - c2)) / 2,
                    ((2 * x - 2) ** 2 * ((c2 + 1) * (x * 2 - 2) + c2) + 2) / 2,
                ),
            },
            "Bounce": {
                "In": 1 - bounce_out(1 - x),
                "Out": bounce_out(x),
                "In Out": np.where(
                    x < 0.5,
                    (1 - bounce_out(1 - 2 * x)) / 2,
                    (1 + bounce_out(2 * x - 1)) / 2,
                ),
            },
            "Elastic": {
                "In": pinned(-(2 ** (10 * x - 10)) * np.sin((x * 10 - 10.75) * c4)),
                "Out": pinned(2 ** (-10 * x) * np.sin((x * 10 - 0.75) * c4) + 1),
                "In Out": pinned(
                    np.where(
                        x < 0.5,
                        -(2 ** (20 * x - 10) * np.sin((20 * x - 11.125) * c5)) / 2,
                        (2 ** (-20 * x + 10) * np.sin((20 * x - 11.125) * c5)) / 2 + 1,
                    )
                ),
            },
        }
    return table[curve][ease]


EASE_CURVES = [
    "Linear",
    "Sinusoidal",
    "Quadratic",
    "Cubic",
    "Quartic",
    "Quintic",
    "Exponential",
    "Circular",
    "Back",
    "Bounce",
    "Elastic",
]
EASE_TYPES = ["In", "Out", "In Out"]


def test_animate_ease_penner():
    mol = mn.Molecule.fetch("4ozs", cache=data_dir)
    grid = np.linspace(0.0, 1.0, 9)
    t = grid[np.arange(len(mol)) % len(grid)]
    mol.store_named_attribute(t, "t")
    combos = [(curve, ease) for curve in EASE_CURVES for ease in EASE_TYPES]

    with mol.tree.reset() as (atoms, join):
        geometry = atoms
        for curve, ease in combos:
            geometry = geometry >> StoreNamedAttribute.point.float(
                name=f"{curve}|{ease}",
                value=(
                    NamedAttribute.float("t").o.attribute
                    >> EaseValue(interpolation=curve, ease=ease)
                ),
            )
        geometry >> join

    for curve, ease in combos:
        got = mol.named_attribute(f"{curve}|{ease}", evaluate=True)
        assert np.allclose(got, _penner(curve, ease, t), atol=1e-4), (curve, ease)

    # Out is the point reflection of In about (0.5, 0.5); the grid is symmetric
    # so both t and 1 - t are sampled.
    for curve in EASE_CURVES:
        ease_in = mol.named_attribute(f"{curve}|In", evaluate=True)
        ease_out = mol.named_attribute(f"{curve}|Out", evaluate=True)
        mirrored = np.array(
            [1 - ease_in[np.argmax(np.isclose(t, 1 - value))] for value in t]
        )
        assert np.allclose(ease_out, mirrored, atol=1e-4), curve


def test_animate_ease_range_and_clamp():
    mol = mn.Molecule.fetch("4ozs", cache=data_dir)
    raw = np.linspace(-1.0, 2.0, 7)[np.arange(len(mol)) % 7]
    mol.store_named_attribute(raw, "raw")

    def ease(clamp: bool):
        return EaseValue(
            value=NamedAttribute.float("raw"),
            interpolation="Linear",
            ease="In",
            clamp=clamp,
            from_=2.0,
            to=-3.0,
        )

    with mol.tree.reset() as (atoms, join):
        (
            atoms
            >> StoreNamedAttribute.point.float(name="clamped", value=ease(True))
            >> StoreNamedAttribute.point.float(name="free", value=ease(False))
            >> join
        )

    clamped = mol.named_attribute("clamped", evaluate=True)
    free = mol.named_attribute("free", evaluate=True)
    assert np.allclose(clamped, 2.0 - 5.0 * np.clip(raw, 0, 1), atol=1e-5)
    assert np.allclose(free, 2.0 - 5.0 * raw, atol=1e-5)


def _stagger(value, rank, span, width):
    # every ID gets an equal window; starts are spread so ~`width` overlap
    denominator = span + width
    start = rank / denominator
    if width == 0:
        return (value >= start).astype(float)
    return np.clip((value - start) / (width / denominator), 0.0, 1.0)


def test_stagger_value():
    mol = mn.Molecule.fetch("4ozs", cache=data_dir)
    index = np.arange(len(mol))
    rank = index % 4
    mol.store_named_attribute(rank, "rank")
    # ids 0..7 split into two groups of 0..3 and 4..7
    mol.store_named_attribute(index % 8, "grouped")
    mol.store_named_attribute(index % 8 // 4, "group")
    values = [0.0, 0.2, 0.5, 0.75, 1.0]
    widths = [0.0, 1.0, 3.0]

    with mol.tree.reset() as (atoms, join):
        geometry = atoms
        for value in values:
            for width in widths:
                geometry = geometry >> StoreNamedAttribute.point.float(
                    name=f"{value}|{width}",
                    value=StaggerValue(
                        value, width=width, id=NamedAttribute.integer("rank")
                    ),
                )
            geometry = geometry >> StoreNamedAttribute.point.float(
                name=f"group|{value}",
                value=StaggerValue(
                    value,
                    width=2.0,
                    id=NamedAttribute.integer("grouped"),
                    group_id=NamedAttribute.integer("group"),
                ),
            )
        geometry >> join

    for value in values:
        for width in widths:
            got = mol.named_attribute(f"{value}|{width}", evaluate=True)
            expected = _stagger(value, rank, 3, width)
            assert np.allclose(got, expected, atol=1e-5), (value, width)
        # each group staggers over its own ID range
        got = mol.named_attribute(f"group|{value}", evaluate=True)
        assert np.allclose(got, _stagger(value, rank, 3, 2.0), atol=1e-5), value


def test_animate_action():
    mol = mn.Molecule.fetch("4ozs", cache=data_dir)
    rank = np.arange(len(mol)) % 4
    mol.store_named_attribute(rank, "rank")
    times = [0.0, 12.0, 17.0, 20.0, 30.0]

    def action(time, **kwargs):
        return AnimateAction(
            start=10.0, length=10.0, time="Value", value=time, **kwargs
        )

    with mol.tree.reset() as (atoms, join):
        geometry = atoms
        for time in times:
            node = action(time)
            geometry = (
                geometry
                >> StoreNamedAttribute.point.float(name=f"linear|{time}", value=node)
                >> StoreNamedAttribute.point.boolean(
                    name=f"active|{time}", value=node.o.active
                )
                >> StoreNamedAttribute.point.float(
                    name=f"cubic|{time}",
                    value=action(time, interpolation="Cubic", ease="In"),
                )
                >> StoreNamedAttribute.point.float(
                    name=f"stagger|{time}",
                    value=action(
                        time,
                        stagger=True,
                        width=1.0,
                        id=NamedAttribute.integer("rank"),
                    ),
                )
                >> StoreNamedAttribute.point.float(
                    name=f"residue|{time}",
                    value=action(time, stagger=True, width=1.0, id=UResID()),
                )
            )
        (
            geometry
            >> StoreNamedAttribute.point.float(name="end", value=action(0.0).o.end)
            >> StoreNamedAttribute.point.float(
                name="frames",
                value=AnimateAction(start=10.0, length=10.0, time="Frames"),
            )
            >> join
        )

    ures_id = mol.named_attribute("ures_id")
    ures_rank = ures_id - ures_id.min()
    for time in times:
        linear = np.clip((time - 10.0) / 10.0, 0.0, 1.0)
        got = mol.named_attribute(f"linear|{time}", evaluate=True)
        assert np.allclose(got, linear, atol=1e-5), time
        active = mol.named_attribute(f"active|{time}", evaluate=True)
        assert np.all(active == (10.0 <= time <= 20.0)), time
        cubic = mol.named_attribute(f"cubic|{time}", evaluate=True)
        assert np.allclose(cubic, linear**3, atol=1e-5), time
        stagger = mol.named_attribute(f"stagger|{time}", evaluate=True)
        assert np.allclose(stagger, _stagger(linear, rank, 3, 1.0), atol=1e-5), time
        residue = mol.named_attribute(f"residue|{time}", evaluate=True)
        expected = _stagger(linear, ures_rank, ures_rank.max(), 1.0)
        assert np.allclose(residue, expected, atol=1e-5), time

    assert np.allclose(mol.named_attribute("end", evaluate=True), 20.0)

    scene = bpy.context.scene
    frame = scene.frame_current
    try:
        scene.frame_set(15)
        frames = mol.named_attribute("frames", evaluate=True)
        assert np.allclose(frames, 0.5, atol=1e-5)
    finally:
        scene.frame_set(frame)


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


def _colormap_menus() -> list[tuple[str, str]]:
    """(category input, colormap) for every map on the Color matplotlib node."""
    import inspect
    from typing import Literal, get_args, get_origin
    from molecularnodes.nodes.geometry import ColorMatplotlib

    params = inspect.signature(ColorMatplotlib.__init__).parameters
    return [
        (param, name)
        for param in ("uniform", "sequential", "sequential_2", "diverging")
        + ("cyclic", "qualitative", "miscellaneous")
        for arg in get_args(params[param].annotation)
        if get_origin(arg) is Literal
        for name in get_args(arg)
    ]


# maps with more detail than 32 stops can hold, in 1/255 of sRGB
_COLORMAP_LOSSY = {
    "gist_ncar": 14.0,
    "nipy_spectral": 11.0,
    "hsv": 7.0,
    "gist_rainbow": 7.0,
    "gist_stern": 4.5,
    "jet": 4.5,
    "gnuplot2": 4.0,
}


def _evaluate_colormap(x: np.ndarray, param: str, name: str, reverse=False):
    """sRGB colours the Color matplotlib node gives at values ``x``."""
    from molecularnodes.nodes.geometry import ColorMatplotlib

    mol = mn.Molecule.fetch("4ozs", cache=data_dir)
    values = np.zeros(len(mol.atoms), dtype=np.float32)
    values[: len(x)] = x
    mol.store_named_attribute(values, "cmap_x")
    with mol.tree.reset() as (atoms, join):
        colormap = ColorMatplotlib(
            value=NamedAttribute.float("cmap_x"),
            reverse=reverse,
            category=param.replace("_", " ").title(),
            **{param: name},
        )
        (
            atoms
            >> StoreNamedAttribute.point.color(name="cmap", value=colormap.o.color)
            >> join
        )
    # the node works in linear RGB, matplotlib in sRGB
    linear = np.clip(mol.named_attribute("cmap", evaluate=True)[: len(x), :3], 0, 1)
    return np.where(
        linear <= 0.0031308, linear * 12.92, 1.055 * linear ** (1 / 2.4) - 0.055
    )


@pytest.mark.parametrize("param, name", _colormap_menus())
def test_color_matplotlib_matches_matplotlib(param, name):
    import matplotlib
    from matplotlib.colors import ListedColormap

    cmap = matplotlib.colormaps[name]
    if isinstance(cmap, ListedColormap) and cmap.N <= 32:
        # constant ramps: sample bin centres, away from the steps
        x = (np.arange(cmap.N) + 0.5) / cmap.N
        tolerance = 0.5
    else:
        # the 256 entries of matplotlib's lookup table
        x = np.linspace(0.0, 1.0, 256)
        tolerance = _COLORMAP_LOSSY.get(name, 2.5)
    srgb = _evaluate_colormap(x, param, name)
    assert np.abs(srgb - cmap(x)[:, :3]).max() * 255 <= tolerance


@pytest.mark.parametrize("mode", ["LINEAR", "B_SPLINE"])
def test_colormap_generator_ramp_model(mode):
    """The generator fits ramps against a numpy model of Blender's Color Ramp;
    check the model against Blender's own evaluation."""
    import importlib.util
    from pathlib import Path

    path = Path(__file__).parents[1] / "docs/dev/colormap_node.py"
    spec = importlib.util.spec_from_file_location("colormap_node", path)
    generator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(generator)

    tree = bpy.data.node_groups.new("ramp_model", "GeometryNodeTree")
    ramp = tree.nodes.new("ShaderNodeValToRGB").color_ramp
    ramp.interpolation = mode
    rng = np.random.default_rng(0)
    x = np.concatenate([np.linspace(0, 1, 101), rng.random(100)])
    for n in range(2, 12):
        pos = np.sort(np.concatenate([[0.0, 1.0], rng.random(n - 2)]))
        colors = rng.random((n, 3))
        while len(ramp.elements) > 1:
            ramp.elements.remove(ramp.elements[-1])
        ramp.elements[0].position = 0.0
        ramp.elements[0].color = (*colors[0], 1.0)
        for p, c in zip(pos[1:], colors[1:]):
            ramp.elements.new(p).color = (*c, 1.0)
        blender = np.array([ramp.evaluate(v)[:3] for v in x])
        model = np.clip(generator.ramp_basis(pos, x, mode) @ colors, 0, 1)
        assert np.abs(blender - model).max() < 1e-3
    bpy.data.node_groups.remove(tree)


def test_color_matplotlib_reverse():
    import matplotlib

    x = np.linspace(0.0, 1.0, 256)
    srgb = _evaluate_colormap(x, "uniform", "viridis", reverse=True)
    expected = matplotlib.colormaps["viridis_r"](x)[:, :3]
    assert np.abs(srgb - expected).max() * 255 <= 2.5


@pytest.mark.parametrize("fade", [0.5, 1.0])
def test_fade_geometry_scales_alpha(fade):
    mol = mn.Molecule.fetch("4ozs", cache=data_dir)
    color = mol.named_attribute("Color")
    with mol.tree.reset() as (atoms, join):
        atoms >> FadeGeometry(fade=fade) >> join

    faded = mol.named_attribute("Color", evaluate=True)
    assert np.allclose(faded[:, :3], color[:, :3])
    assert np.allclose(faded[:, 3], color[:, 3] * fade)


def test_fade_geometry_zero_removes_geometry():
    mol = mn.Molecule.fetch("4ozs", cache=data_dir)
    with mol.tree.reset() as (atoms, join):
        atoms >> FadeGeometry(fade=0.0) >> join

    evaluated = mol.object.evaluated_get(bpy.context.evaluated_depsgraph_get())
    assert len(evaluated.data.vertices) == 0
