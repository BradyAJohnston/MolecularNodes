import random
from typing import Any
import bpy
import MDAnalysis as mda
import numpy as np
import pytest
from databpy.nodes import get_input, get_output
from MDAnalysis.tests.datafiles import DCD, GRO, PSF, XTC
from nodebpy.nodes.geometry import (
    GetBundleItem,
    GetGeometryBundle,
    Group,
    NamedAttribute,
    Points,
    RealizeInstances,
    SetPosition,
    StoreNamedAttribute,
    Value,
)
import molecularnodes as mn
from molecularnodes.nodes._utils import (
    custom_boolean_iswitch,
    custom_color_iswitch,
    get_final_style_nodes,
)
from molecularnodes.nodes.geometry import (
    AnimateEase,
    AnimateReveal,
    AnimateStagger,
    BreakBonds,
    BuildElasticNetwork,
    Charge,
    FindBonds,
    NucleicChi,
    NucleicDihedral,
    PeptideChi,
    PeptideDihedral,
    PeriodicArray,
    SegmentID,
    SetColor,
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
                value=AnimateEase(
                    value=NamedAttribute.float("t"), interpolation=curve, ease=ease
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
        return AnimateEase(
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


def test_animate_stagger():
    mol = mn.Molecule.fetch("4ozs", cache=data_dir)
    rank = (np.arange(len(mol)) % 4).astype(float)
    mol.store_named_attribute(rank, "rank")
    frames = [5.0, 12.0, 20.0, 27.5, 40.0]

    def stagger(frame: float, **kwargs):
        # Frame defaults to the scene frame, so a constant has to be linked in
        options = dict(
            order="Attribute",
            attribute=NamedAttribute.float("rank"),
            frame_start=10,
            delay=5.0,
            length=10.0,
            interpolation="Linear",
            ease="In",
        )
        options.update(kwargs)
        return AnimateStagger(frame=Value(frame), **options)

    with mol.tree.reset() as (atoms, join):
        geometry = atoms
        for frame in frames:
            geometry = geometry >> StoreNamedAttribute.point.float(
                name=f"stagger_{frame}", value=stagger(frame)
            )
        (
            geometry
            >> StoreNamedAttribute.point.float(
                name="start", value=stagger(20.0).o.start
            )
            >> StoreNamedAttribute.point.float(
                name="reverse", value=stagger(20.0, reverse=True)
            )
            >> StoreNamedAttribute.point.float(
                name="snap", value=stagger(20.0, length=0.0)
            )
            >> StoreNamedAttribute.point.float(
                name="cubic", value=stagger(17.0, interpolation="Cubic")
            )
            >> StoreNamedAttribute.point.float(
                name="residue",
                value=AnimateStagger(
                    frame=Value(100.0),
                    frame_start=0,
                    delay=1.0,
                    length=50.0,
                    interpolation="Linear",
                    ease="In",
                ),
            )
            >> join
        )

    start = 10.0 + rank * 5.0
    assert np.allclose(mol.named_attribute("start", evaluate=True), start)
    for frame in frames:
        expected = np.clip((frame - start) / 10.0, 0.0, 1.0)
        got = mol.named_attribute(f"stagger_{frame}", evaluate=True)
        assert np.allclose(got, expected, atol=1e-5), frame

    start_reversed = 10.0 + (rank.max() - rank) * 5.0
    expected = np.clip((20.0 - start_reversed) / 10.0, 0.0, 1.0)
    assert np.allclose(mol.named_attribute("reverse", evaluate=True), expected)

    snap = mol.named_attribute("snap", evaluate=True)
    assert np.array_equal(snap, (20.0 >= start).astype(float))

    linear = np.clip((17.0 - start) / 10.0, 0.0, 1.0)
    cubic = mol.named_attribute("cubic", evaluate=True)
    assert np.allclose(cubic, linear**3, atol=1e-5)

    # the default Order staggers by residue through the `ures_id` attribute
    ures_id = mol.named_attribute("ures_id")
    expected = np.clip((100.0 - ures_id) / 50.0, 0.0, 1.0)
    assert np.allclose(mol.named_attribute("residue", evaluate=True), expected)


def test_animate_reveal():
    mol = mn.Molecule.fetch("4ozs", cache=data_dir)
    n = len(mol)
    factor = (np.arange(n) % 3) / 2.0
    first_half = np.arange(n) < n // 2
    mol.store_named_attribute(factor, "f")
    mol.store_named_attribute(first_half, "first_half")
    color = mol.named_attribute("Color")
    vdw_radii = mol.named_attribute("vdw_radii")
    position = mol.named_attribute("position")
    assert np.all(color[:, 3] == 1.0)

    def reveal(selection=None, **kwargs):
        # nodes must be built inside the tree context, so the selection is a factory
        with mol.tree.reset() as (atoms, join):
            if selection is not None:
                kwargs["selection"] = selection()
            (atoms >> AnimateReveal(factor=NamedAttribute.float("f"), **kwargs) >> join)
        return {
            name: mol.named_attribute(name, evaluate=True)
            for name in ("Color", "vdw_radii", "position")
        }

    alpha = reveal(mode="Alpha")
    assert np.allclose(alpha["Color"][:, :3], color[:, :3])
    assert np.allclose(alpha["Color"][:, 3], factor)
    assert np.allclose(alpha["vdw_radii"], vdw_radii)
    assert np.allclose(alpha["position"], position)

    inverted = reveal(mode="Alpha", invert=True)
    assert np.allclose(inverted["Color"][:, 3], 1.0 - factor)

    # the selection limits the change; unselected points keep their alpha
    selected = reveal(
        mode="Alpha", selection=lambda: NamedAttribute.boolean("first_half")
    )
    assert np.allclose(selected["Color"][first_half, 3], factor[first_half])
    assert np.all(selected["Color"][~first_half, 3] == 1.0)

    # a second reveal multiplies into the alpha the first one wrote
    with mol.tree.reset() as (atoms, join):
        (
            atoms
            >> AnimateReveal(factor=NamedAttribute.float("f"))
            >> AnimateReveal(factor=0.5)
            >> join
        )
    assert np.allclose(mol.named_attribute("Color", evaluate=True)[:, 3], factor / 2)

    scale = reveal(mode="Scale")
    assert np.allclose(scale["vdw_radii"], vdw_radii * factor)
    assert np.allclose(scale["Color"], color)

    cull = reveal(mode="Cull")
    assert len(cull["position"]) == int((factor > 0).sum())
    assert np.allclose(cull["position"], position[factor > 0])
