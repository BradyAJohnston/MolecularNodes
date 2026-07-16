import random
from typing import Any
import bpy
import numpy as np
import pytest
from MDAnalysis.tests.datafiles import DCD, GRO, PSF, XTC
import molecularnodes as mn
from molecularnodes.nodes import nodes
from molecularnodes.nodes.geometry import (
    TopologyBreakBonds,
    TopologyFindBonds,
)
from .constants import codes, data_dir
from .utils import GeometrySet, NumpySnapshotExtension

random.seed(6)


def test_get_nodes():
    mol = mn.Molecule.fetch("4ozs", cache=data_dir).add_style("spheres")

    assert nodes.get_nodes_last_output(mol.node_group)[0].name == "Join Geometry"
    nodes.realize_instances(mol.object)
    assert nodes.get_nodes_last_output(mol.node_group)[0].name == "Realize Instances"
    assert nodes.get_style_node(mol.object).name == "Style Spheres"

    mol2 = mn.Molecule.fetch("1cd3", cache=data_dir).add_style("cartoon", assembly=True)

    assert nodes.get_nodes_last_output(mol2.node_group)[0].name == "Assembly 1cd3"
    assert nodes.get_style_node(mol2.object).name == "Style Cartoon"


def test_selection():
    chain_ids = [let for let in "ABCDEFG123456"]
    node = nodes.custom_boolean_iswitch("test_node", chain_ids, prefix="Chain ")

    input_sockets = nodes.inputs(node)
    for letter, socket in zip(chain_ids, input_sockets.values()):
        assert f"Chain {letter}" == socket.name
        assert socket.default_value is False


@pytest.mark.parametrize("code", codes)
@pytest.mark.parametrize("attribute", ["chain_id", "entity_id"])
def test_selection_working(snapshot_custom: NumpySnapshotExtension, attribute, code):
    mol = mn.Molecule.fetch(code, cache=data_dir).add_style("ribbon")
    group = mol.node_group
    node_sel = nodes.add_selection(
        group, mol.name, mol.object[f"{attribute}s"], attribute
    )

    _n = len(node_sel.inputs)

    nodes.realize_instances(mol.object)

    for inp in node_sel.inputs:
        inp.default_value = True
        pos = mol.named_attribute("position", evaluate=True)
        assert snapshot_custom == pos.shape
        assert snapshot_custom == pos
        inp.default_value = False


@pytest.mark.parametrize("code", codes)
@pytest.mark.parametrize("attribute", ["chain_id", "entity_id"])
def test_color_custom(snapshot_custom: NumpySnapshotExtension, code, attribute):
    mol = mn.Molecule.fetch(code, cache=data_dir).add_style("ribbon")

    group_col = nodes.custom_color_iswitch(
        name=f"Color Entity {mol.name}",
        items=mol.object[f"{attribute}s"],
        attribute_name=attribute,
    )
    group = mol.node_group
    node_col = nodes.add_custom(group, group_col.name, [0, -200])
    group.links.new(node_col.outputs[0], group.nodes["Set Color"].inputs["Color"])

    for i, input in enumerate(node_col.inputs):
        setattr(input, "default_value", mn.color.random_rgb(i))

    assert snapshot_custom == mol.named_attribute("Color")


def test_iswitch_creation():
    items = [str(x) for x in range(10)]
    tree_boolean = nodes.custom_boolean_iswitch("newboolean", items)
    # ensure there isn't an item called 'Color' in the created interface
    assert not tree_boolean.interface.items_tree.get("Color")
    assert tree_boolean.interface.items_tree["Selection"].in_out == "OUTPUT"
    assert tree_boolean.interface.items_tree["Inverted"].in_out == "OUTPUT"
    for i in items:
        assert tree_boolean.interface.items_tree[str(i)].in_out == "INPUT"

    tree_rgba = nodes.custom_color_iswitch("newcolor", items)
    # ensure there isn't an item called 'selection'
    assert not tree_rgba.interface.items_tree.get("Selection")
    assert tree_rgba.interface.items_tree["Color"].in_out == "OUTPUT"
    for i in items:
        assert tree_rgba.interface.items_tree[str(i)].in_out == "INPUT"


def test_op_custom_color():
    mol = mn.Molecule.load(data_dir / "1cd3.cif")
    mol.object.select_set(True)
    group = nodes.custom_color_iswitch(
        name=f"Color Chain {mol.name}", items=mol.object["chain_ids"]
    )

    assert group
    assert group.interface.items_tree["G"].name == "G"
    assert group.interface.items_tree[-1].name == "G"
    assert group.interface.items_tree[0].name == "Color"


def test_color_lookup_supplied():
    col = mn.color.random_rgb(6)
    name = "test"
    node = nodes.custom_color_iswitch(
        name=name,
        items={str(x): col for x in range(10, 20)},
        offset=10,
    )
    assert node.name == name
    for item in nodes.inputs(node).values():
        assert np.allclose(np.array(item.default_value), col)

    node = nodes.custom_color_iswitch(name="test2", items=range(10, 20), offset=10)
    for item in nodes.inputs(node).values():
        assert not np.allclose(np.array(item.default_value), col)


def get_links(sockets):
    for socket in sockets:
        for link in socket.links:
            yield link


@pytest.fixture
def pdb_8h1b():
    return mn.Molecule.fetch("8H1B", cache=data_dir)


# topology_node_names = [n for n in dir(mn.nodes.geometry) if not n.startswith(".")]
topology_node_names = []


@pytest.mark.parametrize("node_name", topology_node_names)
@pytest.mark.parametrize("code", codes)
def test_node_topology(snapshot_custom: NumpySnapshotExtension, code, node_name):
    mol = mn.Molecule.fetch(code, cache=data_dir)

    group = nodes.get_mod(mol.object).node_group = nodes.new_tree()

    group.links.new(
        group.nodes["Group Input"].outputs[0], group.nodes["Group Output"].inputs[0]
    )
    node_att = group.nodes.new("GeometryNodeStoreNamedAttribute")
    node_att.inputs[2].default_value = "test_attribute"
    nodes.insert_last_node(group, node_att)
    # exclude these particular nodes, as they aren't field nodes and so we shouldn't
    # be testing them here. Will create their own particular tests later
    if any(
        keyword in node_name
        for keyword in [
            "Bonds",
            "Bond Count",
            "DSSP",
            "Sample Atomic Attributes",
            "Chain Group ID",
            "Peptide Dihedral",
            "Peptide Chi",
            "Nucleic Dihedral",
            "Nucleic Chi",
        ]
    ):
        return None

    node_topo = nodes.add_custom(
        group, node_name, location=[x - 300 for x in node_att.location]
    )

    if mn.nodes.arrange.node_has_geo_socket(node_topo):
        return None

    if node_name == "Residue Mask":
        node_topo.inputs["atom_name"].default_value = 61

    type_to_data_type = {
        "VECTOR": "FLOAT_VECTOR",
        "VALUE": "FLOAT",
        "BOOLEAN": "BOOLEAN",
        "INT": "INT",
        "RGBA": "FLOAT_COLOR",
        "ROTATION": "QUATERNION",
    }

    for output in node_topo.outputs:
        node_att.data_type = type_to_data_type[output.type]
        input = node_att.inputs["Value"]

        for link in input.links:
            group.links.remove(link)

        group.links.new(output, input)

        assert snapshot_custom == mol.named_attribute("test_attribute", evaluate=True)


@pytest.mark.parametrize(
    "name", ["Peptide Dihedral", "Nucleic Dihedral", "Peptide Chi", "Nucleic Chi"]
)
@pytest.mark.parametrize("code", ["8H1B", "1BNA"])
def test_dihedral_rotations(snapshot_custom: NumpySnapshotExtension, code, name):
    mol = mn.Molecule.fetch(code, cache=data_dir)
    node_sp = mn.nodes.nodes.insert_before(
        mol.modifier_node_tree.nodes["Group Output"], "GeometryNodeSetPosition"
    )
    node = mn.nodes.nodes.insert_before(node_sp.inputs["Position"], name)
    for input in node.inputs:
        if input.name in ["Position", "Selection"]:
            continue
        input.default_value = 1.0
    assert snapshot_custom == mol.named_attribute("position", evaluate=True)[:100]


def test_topo_bonds():
    mol = mn.Molecule.fetch("1BNA", cache=data_dir)
    nodes.get_mod(mol.object).node_group = nodes.new_tree()
    with mol.tree.reset() as (atoms, join):
        atoms >> TopologyBreakBonds(cutoff=0.0) >> join

    # compare the number of edges before and after deleting them with
    bonds = mol.object.data.edges
    no_bonds = mol.evaluate().data.edges
    assert len(bonds) > len(no_bonds)
    assert len(no_bonds) == 0

    # add the node to find the bonds, and ensure the number of bonds pre and post the nodes
    # are the same (other attributes will be different, but for now this is good)
    with mol.tree.reset() as (atoms, join):
        atoms >> TopologyBreakBonds(cutoff=0.0) >> TopologyFindBonds() >> join

    bonds_new = mol.evaluate().data.edges
    assert len(bonds) == len(bonds_new)


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
    assert list(nodes.get_input(tree).outputs.keys()) == ["Atoms", ""]
    assert list(nodes.get_output(tree).inputs.keys()) == ["Geometry", ""]


def test_reuse_node_group():
    mol = mn.Molecule.fetch("4ozs").add_style("spheres")
    tree = bpy.data.node_groups["MN_4ozs"]
    n_nodes = len(tree.nodes)
    bpy.data.objects.remove(mol.object)
    del mol
    assert n_nodes == len(tree.nodes)
    mn.Molecule.fetch("4ozs")
    assert n_nodes == len(tree.nodes)


def _insert_periodic_array(traj: mn.Trajectory):
    node = mn.nodes.nodes.add_custom(traj.modifier_node_tree, "Periodic Array")
    mn.nodes.nodes.insert_last_node(group=traj.tree.tree, node=node)
    return node


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
    traj = mn.Trajectory.load(GRO, XTC, selection="protein")
    node = _insert_periodic_array(traj)

    traj.set_frame(1)
    defaults_0 = _get_node_defaults(node)
    traj.set_frame(10)
    defaults_10 = _get_node_defaults(node)

    dim_idx = slice(1, 7)
    assert not all([x == y for x, y in zip(defaults_0[dim_idx], defaults_10[dim_idx])])
    # the unit cell dimensions ar currently inputs 1..7 for the node as it is setup so
    # we just subset those and check it matches the universe
    assert np.allclose(defaults_10[dim_idx], traj.universe.trajectory.ts.dimensions)

    # for some reason we need to trigger a proper re-evaluation of the GN node tree
    # by saving to a temp file #TODO: look into and try to fix this
    bpy.ops.wm.save_as_mainfile(filepath=str(tmp_path / "example.blend"))
    assert snapshot == GeometrySet(traj.object)


# this topology doesn't have any dimension information so it should just
# update the positions and _attempt_ to update the periodic box but fail not do so quietly
# and everything remains 0
def test_periodic_array_no_dimensions():
    traj = mn.Trajectory.load(PSF, DCD, selection="protein")
    node = _insert_periodic_array(traj)

    traj.set_frame(1)
    defaults_0 = _get_node_defaults(node)
    traj.set_frame(frame=10)
    defaults_10 = _get_node_defaults(node)

    assert defaults_0 == defaults_10
    assert defaults_0[1:7] == [0] * 6
