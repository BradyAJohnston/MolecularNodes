# import bpy
import numpy as np
import pytest
import molecularnodes as mn
from molecularnodes.blender import nodes
import random

from .utils import sample_attribute, NumpySnapshotExtension
from .constants import codes, data_dir

random.seed(6)

mn._test_register()


def test_node_name_format():
    assert mn.blender.nodes.format_node_name("Style Cartoon") == "Style Cartoon"
    assert (
        mn.blender.nodes.format_node_name("MN_dna_double_helix") == "DNA Double Helix"
    )
    assert (
        mn.blender.nodes.format_node_name("MN_topo_vector_angle")
        == "Topology Vector Angle"
    )


def test_get_nodes():
    bob = mn.io.fetch("4ozs", style="spheres", cache_dir=data_dir).object

    assert (
        nodes.get_nodes_last_output(bob.modifiers["MolecularNodes"].node_group)[0].name
        == "Style Spheres"
    )
    nodes.realize_instances(bob)
    assert (
        nodes.get_nodes_last_output(bob.modifiers["MolecularNodes"].node_group)[0].name
        == "Realize Instances"
    )
    assert nodes.get_style_node(bob).name == "Style Spheres"

    bob2 = mn.io.fetch(
        "1cd3", style="cartoon", build_assembly=True, cache_dir=data_dir
    ).object

    assert (
        nodes.get_nodes_last_output(bob2.modifiers["MolecularNodes"].node_group)[0].name
        == "Assembly 1cd3"
    )
    assert nodes.get_style_node(bob2).name == "Style Cartoon"


def test_selection():
    chain_ids = [let for let in "ABCDEFG123456"]
    node = nodes.custom_iswitch(
        "test_node", chain_ids, prefix="Chain ", dtype="BOOLEAN"
    )

    input_sockets = nodes.inputs(node)
    for letter, socket in zip(chain_ids, input_sockets.values()):
        assert f"Chain {letter}" == socket.name
        assert socket.default_value is False


@pytest.mark.parametrize("code", codes)
@pytest.mark.parametrize("attribute", ["chain_id", "entity_id"])
def test_selection_working(snapshot_custom: NumpySnapshotExtension, attribute, code):
    mol = mn.io.fetch(code, style="ribbon", cache_dir=data_dir).object
    group = mol.modifiers["MolecularNodes"].node_group
    node_sel = nodes.add_selection(group, mol.name, mol[f"{attribute}s"], attribute)

    n = len(node_sel.inputs)

    for i in random.sample(list(range(n)), max(n - 2, 1)):
        node_sel.inputs[i].default_value = True

    nodes.realize_instances(mol)

    assert snapshot_custom == sample_attribute(mol, "position", evaluate=True)


@pytest.mark.parametrize("code", codes)
@pytest.mark.parametrize("attribute", ["chain_id", "entity_id"])
def test_color_custom(snapshot_custom: NumpySnapshotExtension, code, attribute):
    mol = mn.io.fetch(code, style="ribbon", cache_dir=data_dir)

    group_col = mn.blender.nodes.custom_iswitch(
        name=f"Color Entity {mol.name}",
        iter_list=mol.object[f"{attribute}s"],
        field=attribute,
        dtype="RGBA",
    )
    group = mol.object.modifiers["MolecularNodes"].node_group
    node_col = mn.blender.nodes.add_custom(group, group_col.name, [0, -200])
    group.links.new(node_col.outputs[0], group.nodes["Set Color"].inputs["Color"])
    for i, input in enumerate(node_col.inputs):
        input.default_value = mn.color.random_rgb(i)

    assert snapshot_custom == mol.get_attribute("Color")


def test_custom_resid_selection():
    node = mn.blender.nodes.resid_multiple_selection("new_node", "1, 5, 10-20, 40-100")
    numbers = [1, 5, 10, 20, 40, 100]
    assert len(nodes.outputs(node)) == 2
    counter = 0
    for item in node.interface.items_tree:
        if item.in_out == "INPUT":
            assert item.default_value == numbers[counter]
            counter += 1


def test_iswitch_creation():
    items = list(range(10))
    tree_boolean = mn.blender.nodes.custom_iswitch(
        name="newboolean", iter_list=items, dtype="BOOLEAN"
    )
    # ensure there isn't an item called 'Color' in the created interface
    assert not tree_boolean.interface.items_tree.get("Color")
    assert tree_boolean.interface.items_tree["Selection"].in_out == "OUTPUT"
    assert tree_boolean.interface.items_tree["Inverted"].in_out == "OUTPUT"
    for i in items:
        assert tree_boolean.interface.items_tree[str(i)].in_out == "INPUT"

    tree_rgba = mn.blender.nodes.custom_iswitch(
        name="newcolor", iter_list=items, dtype="RGBA"
    )
    # ensure there isn't an item called 'selection'
    assert not tree_rgba.interface.items_tree.get("Selection")
    assert tree_rgba.interface.items_tree["Color"].in_out == "OUTPUT"
    for i in items:
        assert tree_rgba.interface.items_tree[str(i)].in_out == "INPUT"


def test_op_custom_color():
    mol = mn.io.load_local(data_dir / "1cd3.cif").object
    mol.select_set(True)
    group = mn.blender.nodes.custom_iswitch(
        name=f"Color Chain {mol.name}", iter_list=mol["chain_ids"], dtype="RGBA"
    )

    assert group
    assert group.interface.items_tree["G"].name == "G"
    assert group.interface.items_tree[-1].name == "G"
    assert group.interface.items_tree[0].name == "Color"


def test_color_lookup_supplied():
    col = mn.color.random_rgb(6)
    name = "test"
    node = mn.blender.nodes.custom_iswitch(
        name=name,
        iter_list=range(10, 20),
        dtype="RGBA",
        default_values=[col for i in range(10)],
        start=10,
    )
    assert node.name == name
    for item in nodes.inputs(node).values():
        assert np.allclose(np.array(item.default_value), col)

    node = mn.blender.nodes.custom_iswitch(
        name="test2", iter_list=range(10, 20), dtype="RGBA", start=10
    )
    for item in nodes.inputs(node).values():
        assert not np.allclose(np.array(item.default_value), col)


def get_links(sockets):
    links = []
    for socket in sockets:
        for link in socket.links:
            yield link


def test_change_style():
    model = mn.io.fetch("1cd3", style="cartoon", cache_dir=data_dir).object
    style_node_1 = nodes.get_style_node(model).name
    mn.blender.nodes.change_style_node(model, "ribbon")
    style_node_2 = nodes.get_style_node(model).name

    assert style_node_1 != style_node_2

    styles_to_check = ["ribbon", "cartoon", "ball_and_stick", "surface"] + list(
        [f"preset_{i}" for i in [1, 2, 3, 4]]
    )

    for style in styles_to_check:
        style_node_1 = nodes.get_style_node(model)
        links_in_1 = [link.from_socket.name for link in get_links(style_node_1.inputs)]
        links_out_1 = [
            link.from_socket.name for link in get_links(style_node_1.outputs)
        ]

        nodes.change_style_node(model, style)
        style_node_2 = nodes.get_style_node(model)
        links_in_2 = [link.from_socket.name for link in get_links(style_node_2.inputs)]
        links_out_2 = [
            link.from_socket.name for link in get_links(style_node_2.outputs)
        ]

        assert len(links_in_1) == len(links_in_2)
        assert len(links_out_1) == len(links_out_2)


@pytest.fixture
def pdb_8h1b():
    return mn.io.fetch("1bna", del_solvent=False, cache_dir=data_dir)


node_names = [
    node["name"]
    for node in mn.ui.node_info.menu_items["topology"]
    if not node == "break"
]


def test_nodes_exist():
    for item in mn.ui.node_info.menu_items:
        if isinstance(item, str):
            continue
        if hasattr(item, "function"):
            continue
        for name in [node.name for node in item.items()]:
            mn.blender.nodes.append(name)
            assert True


@pytest.mark.parametrize("node_name", node_names)
def test_node_topology(snapshot_custom: NumpySnapshotExtension, pdb_8h1b, node_name):
    mol = pdb_8h1b.object

    group = nodes.get_mod(mol).node_group

    group.links.new(
        group.nodes["Group Input"].outputs[0], group.nodes["Group Output"].inputs[0]
    )
    node_att = group.nodes.new("GeometryNodeStoreNamedAttribute")
    node_att.inputs[2].default_value = "test_attribute"
    nodes.insert_last_node(group, node_att)
    # exclude these particular nodes, as they aren't field nodes and so we shouldn't
    # be testing them here. Will create their own particular tests later
    if any(
        keyword in node_name for keyword in ["Backbone", "Bonds", "Bond Count", "DSSP"]
    ):
        return None

    node_topo = nodes.add_custom(
        group, node_name, location=[x - 300 for x in node_att.location]
    )

    if node_name == "Point Group Mask":
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

        assert snapshot_custom == mn.blender.obj.get_attribute(
            mol, "test_attribute", evaluate=True
        )


def test_compute_backbone(snapshot_custom: NumpySnapshotExtension):
    mol = mn.io.fetch("1CCN", del_solvent=False, cache_dir=data_dir).object

    group = nodes.get_mod(mol).node_group

    group.links.new(
        group.nodes["Group Input"].outputs[0], group.nodes["Group Output"].inputs[0]
    )
    node_att = group.nodes.new("GeometryNodeStoreNamedAttribute")
    node_att.inputs[2].default_value = "test_attribute"
    node_backbone = nodes.add_custom(group, "Topology Compute Backbone")
    nodes.insert_last_node(group, node_backbone)
    nodes.insert_last_node(group, node_att)
    node_names = ["Backbone Positions"]
    for node_name in node_names:
        node_topo = nodes.add_custom(
            group, node_name, location=[x - 300 for x in node_att.location]
        )

        if node_name == "Point Group Mask":
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

            assert snapshot_custom == mn.blender.obj.get_attribute(
                mol, "test_attribute", evaluate=True
            )

        for angle in ["Phi", "Psi"]:
            output = node_backbone.outputs[angle]
            node_att.data_type = type_to_data_type[output.type]
            input = node_att.inputs["Value"]

            for link in input.links:
                group.links.remove(link)

            group.links.new(output, input)

            assert snapshot_custom == mn.blender.obj.get_attribute(
                mol, "test_attribute", evaluate=True
            )


def test_topo_bonds():
    mol = mn.io.fetch("1BNA", del_solvent=True, style=None, cache_dir=data_dir).object
    group = nodes.get_mod(mol).node_group = nodes.new_group()

    # add the node that will break bonds, set the cutoff to 0
    node_break = nodes.add_custom(group, "Topology Break Bonds")
    nodes.insert_last_node(group, node=node_break)
    node_break.inputs["Cutoff"].default_value = 0

    # compare the number of edges before and after deleting them with
    bonds = mol.data.edges
    no_bonds = mn.blender.obj.evaluated(mol).data.edges
    assert len(bonds) > len(no_bonds)
    assert len(no_bonds) == 0

    # add the node to find the bonds, and ensure the number of bonds pre and post the nodes
    # are the same (other attributes will be different, but for now this is good)
    node_find = nodes.add_custom(group, "Topology Find Bonds")
    nodes.insert_last_node(group, node=node_find)
    bonds_new = mn.blender.obj.evaluated(mol).data.edges
    assert len(bonds) == len(bonds_new)
