import bpy
import pytest
import molecularnodes as mn
from .constants import data_dir


def test_add_style_with_selection():
    mol = mn.Molecule.fetch("4ozs").add_style("cartoon")
    mol.store_named_attribute(mol.named_attribute("is_side_chain"), "show_side_chains")
    mol.add_style("ball_and_stick", selection="show_side_chains")

    node_style = [
        node
        for node in mol.modifier_node_tree.nodes
        if (
            isinstance(node, bpy.types.GeometryNodeGroup)
            and node.node_tree.name == "Style Ball and Stick"
        )
    ][0]

    assert (
        node_style.inputs["Selection"].links[0].from_node.inputs["Name"].default_value
        == "show_side_chains"
    )

    with pytest.warns(UserWarning):
        mol.add_style("cartoon", selection="non_existing_selection")


def test_styles_panel_style_node_selection():
    """The Styles panel identifies the active style node from the list index."""
    from molecularnodes.ui.panel import is_style_node

    mol = mn.Molecule.load(data_dir / "1cd3.cif").add_style("cartoon")
    node_group = mol.modifier_node_tree
    style_names = [n.name for n in node_group.nodes if is_style_node(n)]
    assert style_names, "adding a style should create a style node"

    # the panel treats styles_active_index as an index into node_group.nodes and
    # only shows details when it points at a style node
    index = node_group.nodes.find(style_names[0])
    mol.object.mn.styles_active_index = index

    active = node_group.nodes[mol.object.mn.styles_active_index]
    assert is_style_node(active)
    assert active.name == style_names[0]


def test_styles_panel_node_group_lookup_without_session():
    """Style nodes are found from the object alone, without a session link."""
    from molecularnodes.ui.panel import get_entity_node_group, is_style_node

    session = mn.session.get_session()
    mol = mn.Molecule.load(data_dir / "1cd3.cif").add_style("cartoon")
    obj = mol.object

    # drop the session link (as after re-opening a .blend without the pickle)
    session.remove_entity(obj.uuid)
    assert session.get(obj.uuid) is None

    # the node group and its style nodes are still reachable from the object
    node_group = get_entity_node_group(obj)
    assert node_group is not None
    style_names = [n.name for n in node_group.nodes if is_style_node(n)]
    assert style_names


def test_swap_style_operator():
    """mn.swap_style swaps the style node's tree in place, keeping connections."""
    from molecularnodes.ui.panel import is_style_node

    mol = mn.Molecule.load(data_dir / "1cd3.cif").add_style("cartoon")
    ng = mol.modifier_node_tree
    style = [n for n in ng.nodes if is_style_node(n)][0]
    assert style.node_tree.name == "Style Cartoon"

    res = bpy.ops.mn.swap_style(
        name_tree=ng.name, name_node=style.name, style="surface"
    )
    assert res == {"FINISHED"}

    # still exactly one style node, now referencing the Surface style tree
    style_nodes = [n for n in ng.nodes if is_style_node(n)]
    assert len(style_nodes) == 1
    assert style_nodes[0].node_tree.name == "Style Surface"
