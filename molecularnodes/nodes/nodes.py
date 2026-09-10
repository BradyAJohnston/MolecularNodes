from typing import Iterable, Literal
import bpy
from bpy.types import GeometryNodeTree
from databpy.nodes import (
    append_from_blend,
    swap_tree,
)
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from .. import color
from ..assets import MN_DATA_FILE
from . import geometry as mng

NODE_WIDTH = 180

NODE_SPACING = 250


socket_types = {
    "BOOLEAN": "NodeSocketBool",
    "GEOMETRY": "NodeSocketGeometry",
    "INT": "NodeSocketInt",
    "MATERIAL": "NodeSocketMaterial",
    "VECTOR": "NodeSocketVector",
    "STRING": "NodeSocketString",
    "VALUE": "NodeSocketFloat",
    "COLLECTION": "NodeSocketCollection",
    "TEXTURE": "NodeSocketTexture",
    "COLOR": "NodeSocketColor",
    "RGBA": "NodeSocketColor",
    "IMAGE": "NodeSocketImage",
}


STYLE_NODE_MAPPING = {
    "spheres": mng.StyleSpheres,
    "cartoon": mng.StyleCartoon,
    "ribbon": mng.StyleRibbon,
    "surface": mng.StyleSurface,
    "sticks": mng.StyleSticks,
    "ball_and_stick": mng.StyleBallAndStick,
}

STYLE_LITERALS = Literal[
    "spheres", "cartoon", "ribbon", "surface", "sticks", "ball_and_stick"
]

# current implemented representations
styles_mapping = {
    "preset_1": "Style Preset 1",
    "preset_2": "Style Preset 2",
    "preset_3": "Style Preset 3",
    "preset_4": "Style Preset 4",
    "atoms": "Style Spheres",
    "spheres": "Style Spheres",
    "vdw": "Style Spheres",
    "sphere": "Style Spheres",
    "cartoon": "Style Cartoon",
    "sticks": "Style Sticks",
    "ribbon": "Style Ribbon",
    "surface": "Style Surface",
    "ball_and_stick": "Style Ball and Stick",
    "ball+stick": "Style Ball and Stick",
    "oxdna": "MN_oxdna_style_ribbon",
    "density_surface": "Density Style Surface",
    "density_iso_surface": "Density Style ISO Surface",
    "density_wire": "Density Style Wire",
}


def inputs(node):
    items = {}
    for item in node.interface.items_tree:
        if item.item_type == "SOCKET":
            if item.in_out == "INPUT":
                items[item.name] = item
    return items


def create_debug_group(name="MolecularNodesDebugGroup"):
    group = new_tree(name=name, fallback=False)
    info = group.nodes.new("GeometryNodeObjectInfo")
    group.links.new(info.outputs["Geometry"], group.nodes["Group Output"].inputs[0])
    return group


def get_output(group) -> bpy.types.GeometryNode:
    return group.nodes[
        bpy.app.translations.pgettext_data(
            "Group Output",
        )
    ]


def get_input(group) -> bpy.types.GeometryNode:
    return group.nodes[
        bpy.app.translations.pgettext_data(
            "Group Input",
        )
    ]


def get_mod(object, name="Molecular Nodes"):
    node_mod = object.modifiers.get(name)
    if not node_mod:
        node_mod = object.modifiers.new(name, "NODES")
    object.modifiers.active = node_mod
    return node_mod


def get_nodes_last_output(group):
    output = get_output(group)
    last = output.inputs[0].links[0].from_node
    return last, output


def previous_node(node):
    "Get the node which is the first connection to the first input of this node"
    prev = node.inputs[0].links[0].from_node
    return prev


def node_group_name(node) -> str:
    """
    The name of the node group a node instances, or "" if it doesn't instance one.

    The node's own name is not a reliable identifier: swapping the node group of an
    existing node leaves the old name behind, and nodes created via the node API are
    named generically ("Group"). The node group name always tracks what is being used.
    """
    tree = getattr(node, "node_tree", None)
    return tree.name if tree is not None else ""


def style_node(group):
    prev = previous_node(get_output(group))
    while "Style" not in node_group_name(prev):
        prev = previous_node(prev)
    return prev


def get_style_node(object):
    "Walk back through the primary node connections until you find the first style node"
    group = object.modifiers["Molecular Nodes"].node_group
    return style_node(group)


def star_node(group):
    prev = previous_node(get_output(group))
    while "Starfile Instances" not in node_group_name(prev):
        prev = previous_node(prev)
    return prev


def get_star_node(object):
    "Walk back through the primary node connections until you find the first style node"
    group = object.modifiers["Molecular Nodes"].node_group
    return star_node(group)


def insert_last_node(group, node, link_input=True):
    last, output = get_nodes_last_output(group)
    link = group.links.new
    location = output.location
    output.location = [location[0] + 300, location[1]]
    node.location = [location[0] - 300, location[1]]
    if link_input:
        link(last.outputs[0], node.inputs[0])
    link(node.outputs[0], output.inputs[0])


def realize_instances(obj):
    group = obj.modifiers["Molecular Nodes"].node_group
    realize = group.nodes.new("GeometryNodeRealizeInstances")
    insert_last_node(group, realize)


def append(name: str, link: bool = False) -> bpy.types.GeometryNodeTree:
    "Append a GN node from the MN data file"
    GN_TREES_PATH = MN_DATA_FILE / "NodeTree"
    return append_from_blend(name, filepath=str(GN_TREES_PATH), link=link)


def swap(node: bpy.types.Node, tree: str | bpy.types.NodeTree) -> None:
    "Swap out the node's node_tree, maintaining the old socket connections"
    if isinstance(tree, str):
        try:
            tree = bpy.data.node_groups[tree]
        except KeyError:
            tree = append(tree)
    # only change the label if it hasn't been customised away from the tree name
    if node.label == node.node_tree.name:
        node.label = tree.name
    swap_tree(node=node, tree=tree)


def micrograph_material():
    """
    Append MN_micrograph_material to the .blend file it it doesn't already exist,
    and return that material.
    """

    mat_name = "MN_micrograph_material"

    return bpy.data.materials[mat_name]


def new_tree(
    name: str = "Geometry Nodes",
    geometry: bool = True,
    input_name: str = "Geometry",
    output_name: str = "Geometry",
    is_modifier: bool = False,
    fallback: bool = True,
) -> bpy.types.GeometryNodeTree:
    tree = bpy.data.node_groups.get(name)
    # if the group already exists, return it and don't create a new one
    if tree and fallback:
        if not isinstance(tree, bpy.types.GeometryNodeTree):
            raise TypeError(f"Expected a GeometryNodeTree, got {type(tree)}")
        return tree

    # create a new group for this particular name and do some initial setup
    tree: bpy.types.GeometryNodeTree = bpy.data.node_groups.new(
        name=name,
        type="GeometryNodeTree",
    )
    input_node = tree.nodes.new("NodeGroupInput")
    output_node = tree.nodes.new("NodeGroupOutput")
    input_node.location.x = -200 - input_node.width
    output_node.location.x = 200
    if geometry:
        tree.interface.new_socket(
            input_name, in_out="INPUT", socket_type="NodeSocketGeometry"
        )
        tree.interface.new_socket(
            output_name, in_out="OUTPUT", socket_type="NodeSocketGeometry"
        )
        tree.links.new(output_node.inputs[0], input_node.outputs[0])
    tree.is_modifier = is_modifier
    return tree


def last_node(tree: bpy.types.GeometryNodeTree) -> bpy.types.GeometryNode:
    output = get_output(tree)
    try:
        return output.inputs[0].links[0].from_socket.node
    except IndexError:
        return output


def node_previous(node):
    return node.inputs[0].links[0].from_socket.node


def final_join(tree: bpy.types.GeometryNodeTree) -> bpy.types.GeometryNode:
    """
    Get the last JoinGeometry node in the tree.
    """
    # output = get_output(tree)
    current = last_node(tree)
    try:
        while True:
            if current.bl_idname == "GeometryNodeGroupInput":
                raise RuntimeError
            if current.bl_idname == "GeometryNodeJoinGeometry":
                return current
            current = node_previous(current)
    except (RuntimeError, IndexError):
        pass


def custom_boolean_iswitch(
    name: str,
    items: Iterable[str],
    attribute_name: str = "chain_id",
    offset: int = 0,
    prefix: str = "",
) -> TreeBuilder[GeometryNodeTree]:
    with g.tree(name) as tree:
        attr = g.NamedAttribute.integer(attribute_name)

        switch = g.IndexSwitch.boolean(
            index=attr if offset == 0 else attr + offset,
            items=[tree.inputs.boolean(prefix + x) for x in items],
        )

        switch >> tree.outputs.boolean("Selection")
        ~switch >> tree.outputs.boolean("Inverted")

    tree.tree.color_tag = "INPUT"

    return tree


def custom_color_iswitch(
    name: str,
    items: dict[str, tuple[float, float, float, float]] | Iterable[int | float | str],
    attribute_name: str = "chain_id",
    offset: int = 0,
) -> TreeBuilder[GeometryNodeTree]:
    with g.tree(name) as tree:
        attr = g.NamedAttribute.integer(attribute_name)

        if not isinstance(items, dict):
            items: dict[str, tuple[float, ...]] = {
                str(key): color.random_rgb() for key in items
            }

        switch = g.IndexSwitch.color(
            index=attr if offset == 0 else attr + offset,
            items=[tree.inputs.color(key, value) for key, value in items.items()],
        )

        switch >> tree.outputs.color()

    tree.tree.color_tag = "INPUT"
    return tree
