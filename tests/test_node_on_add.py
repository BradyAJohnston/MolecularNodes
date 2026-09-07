"""
Tests for the on-add setup that runs when MN node group assets are imported.

Asset drag-and-drop can't be simulated headlessly, but the drop operator's
steps can: it appends the node group (which fires ``blend_import_post``) and
then creates the node instance in the edited tree. The deferred node setup
normally runs via a timer, which doesn't fire in background mode, so tests
invoke ``_process_pending`` directly.
"""

import bpy
import molecularnodes as mn
from molecularnodes.assets import MN_DATA_FILE
from molecularnodes.nodes import handlers as node_handlers
from .constants import data_dir


def _append_asset(name: str, reuse: bool = False) -> bpy.types.GeometryNodeTree:
    "Append a node group the way an asset drop does, firing blend_import_post."
    bpy.ops.wm.append(
        filepath=f"{MN_DATA_FILE}/NodeTree/{name}",
        directory=f"{MN_DATA_FILE}/NodeTree/",
        filename=name,
        do_reuse_local_id=reuse,
    )
    return bpy.data.node_groups[name]


def _material_socket(tree):
    return next(
        (
            s
            for s in tree.interface.items_tree
            if s.item_type == "SOCKET"
            and s.in_out == "INPUT"
            and s.socket_type == "NodeSocketMaterial"
        ),
        None,
    )


def test_handler_registered():
    assert node_handlers.node_asset_import_post in bpy.app.handlers.blend_import_post


def test_style_import_sets_material_default():
    tree = _append_asset("Style Spheres")
    socket = _material_socket(tree)
    assert socket.default_value is not None
    assert socket.default_value.name == "MN Default"

    # a node instance created from the tree inherits the default, which is
    # what the drop operator does right after the append
    host = bpy.data.node_groups.new("host", "GeometryNodeTree")
    node = host.nodes.new("GeometryNodeGroup")
    node.node_tree = tree
    assert node.inputs["Material"].default_value.name == "MN Default"


def test_style_import_keeps_existing_material_default():
    mat = bpy.data.materials.new("Existing")
    tree = _append_asset("Style Ribbon")
    socket = _material_socket(tree)
    socket.default_value = mat
    # a re-drop that reuses the local group must not clobber the user's choice
    _append_asset("Style Ribbon", reuse=True)
    assert socket.default_value == mat


def test_assembly_instance_drop_builds_data_object():
    node_handlers._pending.clear()
    mol = mn.Molecule.load(data_dir / "1cd3.cif").add_style("cartoon")
    assert mol.assemblies(as_array=True) is not None

    tree = _append_asset("Assembly Instance")
    assert node_handlers._pending

    # the drop operator creates the node after the append
    node = mol.modifier_node_tree.nodes.new("GeometryNodeGroup")
    node.node_tree = tree
    assert node.inputs["Data Object"].default_value is None

    node_handlers._process_pending()
    data_obj = node.inputs["Data Object"].default_value
    assert data_obj is not None
    assert data_obj.name == f".data_{mol.name}_assemblies"
    assert node.get(node_handlers._MARKER)

    # a second drop reuses the same data object and doesn't touch the first node
    _append_asset("Assembly Instance", reuse=True)
    second = mol.modifier_node_tree.nodes.new("GeometryNodeGroup")
    second.node_tree = tree
    node_handlers._process_pending()
    assert second.inputs["Data Object"].default_value == data_obj
    assert node.inputs["Data Object"].default_value == data_obj


def test_assembly_instance_without_entity_is_untouched():
    node_handlers._pending.clear()
    # dropped into a tree that belongs to no MN entity: setup must no-op
    tree = _append_asset("Assembly Instance")
    host = bpy.data.node_groups.new("host", "GeometryNodeTree")
    node = host.nodes.new("GeometryNodeGroup")
    node.node_tree = tree

    node_handlers._process_pending()
    assert node.inputs["Data Object"].default_value is None
    assert not node_handlers._pending
