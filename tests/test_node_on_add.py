"""
Tests for the on-add setup that runs when MN node group assets are imported.

Asset drag-and-drop can't be simulated headlessly, but the add operator's
steps can: it imports the node group (which fires ``blend_import_post``) and
then creates the node instance in the edited tree. Depending on the asset
library's import method the group is appended or linked (the default "Pack"
method is a link plus packing), so both paths are covered. The deferred node
setup normally runs via a timer, which doesn't fire in background mode, so
tests invoke ``_process_pending`` directly.
"""

import bpy
import pytest
import molecularnodes as mn
from molecularnodes.assets import MN_DATA_FILE
from molecularnodes.nodes import handlers as node_handlers
from .constants import data_dir


def _import_asset(
    name: str, link: bool = False, reuse: bool = False
) -> bpy.types.GeometryNodeTree:
    "Import a node group the way an asset add does, firing blend_import_post."
    op = bpy.ops.wm.link if link else bpy.ops.wm.append
    kwargs = {} if link else {"do_reuse_local_id": reuse}
    op(
        filepath=f"{MN_DATA_FILE}/NodeTree/{name}",
        directory=f"{MN_DATA_FILE}/NodeTree/",
        filename=name,
        **kwargs,
    )
    return next(
        ng
        for ng in bpy.data.node_groups
        if ng.name == name and (ng.library is not None) == link
    )


def _new_group_node(tree, host=None):
    if host is None:
        host = bpy.data.node_groups.new("host", "GeometryNodeTree")
    node = host.nodes.new("GeometryNodeGroup")
    node.node_tree = tree
    return node


def test_handler_registered():
    assert node_handlers.node_asset_import_post in bpy.app.handlers.blend_import_post


@pytest.mark.parametrize("link", [False, True], ids=["appended", "linked"])
def test_style_import_sets_node_material(link):
    node_handlers._pending.clear()
    tree = _import_asset("Style Spheres", link=link)
    assert node_handlers._pending

    node = _new_group_node(tree)
    assert node.inputs["Material"].default_value is None

    node_handlers._process_pending()
    assert node.inputs["Material"].default_value is not None
    assert node.inputs["Material"].default_value.name == "MN Default"
    assert node.get(node_handlers._MARKER)


def test_style_import_keeps_existing_material():
    node_handlers._pending.clear()
    mat = bpy.data.materials.new("Existing")
    tree = _import_asset("Style Ribbon")
    node = _new_group_node(tree)
    node.inputs["Material"].default_value = mat

    # a re-add that reuses the group must not clobber the already-set value
    _import_asset("Style Ribbon", reuse=True)
    node_handlers._process_pending()
    assert node.inputs["Material"].default_value == mat


@pytest.mark.parametrize("link", [False, True], ids=["appended", "linked"])
def test_assembly_instance_add_builds_data_object(link):
    node_handlers._pending.clear()
    mol = mn.Molecule.load(data_dir / "1cd3.cif").add_style("cartoon")
    assert mol.assemblies(as_array=True) is not None

    tree = _import_asset("Assembly Instance", link=link)
    assert node_handlers._pending

    # the add operator creates the node after the import
    node = _new_group_node(tree, host=mol.modifier_node_tree)
    assert node.inputs["Data Object"].default_value is None

    node_handlers._process_pending()
    data_obj = node.inputs["Data Object"].default_value
    assert data_obj is not None
    assert data_obj.name == f".data_{mol.name}_assemblies"
    assert node.get(node_handlers._MARKER)

    # a second add reuses the same data object and doesn't touch the first node
    _import_asset("Assembly Instance", link=link, reuse=True)
    second = _new_group_node(tree, host=mol.modifier_node_tree)
    node_handlers._process_pending()
    assert second.inputs["Data Object"].default_value == data_obj
    assert node.inputs["Data Object"].default_value == data_obj


def test_assembly_instance_without_entity_is_untouched():
    node_handlers._pending.clear()
    # added to a tree that belongs to no MN entity: setup must no-op
    tree = _import_asset("Assembly Instance")
    node = _new_group_node(tree)

    node_handlers._process_pending()
    assert node.inputs["Data Object"].default_value is None
    assert not node_handlers._pending
