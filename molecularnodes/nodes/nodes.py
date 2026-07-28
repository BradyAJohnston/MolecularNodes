from typing import Iterable, Literal
import bpy
from databpy.nodes import (
    append_from_blend,
)
from mathutils import Vector
from nodebpy import geometry as g
from .. import color, utils
from ..assets import MN_DATA_FILE
from ..blender import mesh
from . import geometry as mng
from .material import assign_material

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


def set_selection(group, node, selection):
    pos = node.location
    pos = [pos[0] - 200, pos[1] - 200]
    selection.location = pos
    group.links.new(selection.outputs[0], node.inputs["Selection"])

    return selection


def create_debug_group(name="MolecularNodesDebugGroup"):
    group = new_tree(name=name, fallback=False)
    info = group.nodes.new("GeometryNodeObjectInfo")
    group.links.new(info.outputs["Geometry"], group.nodes["Group Output"].inputs[0])
    return group


def add_selection(group, sel_name, input_list, field="chain_id"):
    style = style_node(group)
    sel_node = add_custom(
        group,
        custom_boolean_iswitch(
            name="selection", items=input_list, attribute_name=field
        ).name,
    )

    set_selection(group, style, sel_node)
    return sel_node


def get_selection(node: bpy.types.GeometryNode) -> bpy.types.GeometryNode | None:
    sel_input = node.inputs.get("Selection")
    if not sel_input:
        return None
    try:
        return sel_input.links[0].from_socket.node
    except (KeyError, IndexError):
        return None


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


def add_custom(
    group: bpy.types.GeometryNodeTree,
    name: str,
    location: list[float, float] | Vector = [0, 0],
    width: float = NODE_WIDTH,
    material: str | bpy.types.Material = "default",
    show_options: bool = False,
    link: bool = False,
) -> bpy.types.GeometryNodeGroup:
    node: bpy.types.GeometryNodeGroup = group.nodes.new("GeometryNodeGroup")
    node.node_tree = append(name, link=link)
    # set the label to the node tree name by default
    node.label = node.node_tree.name

    # if there is an input socket called 'Material', assign it to the base MN material
    # if another material is supplied, use that instead.
    assign_material(node, new_material=material)

    # move and format the node for arranging
    node.location = location
    node.width = width
    node.show_options = show_options
    node.name = name

    return node


def assembly_data_object_from_obj(obj: bpy.types.Object) -> bpy.types.Object:
    data_obj_name = f".data_{obj.name}_assemblies"
    data_obj = bpy.data.objects.get(data_obj_name)
    if not data_obj:
        transforms = utils.array_quaternions_from_dict(obj.mn.biological_assemblies)
        data_obj = mesh.create_data_object(array=transforms, name=data_obj_name)

    return data_obj


def insert_join_last(tree: bpy.types.GeometryNodeTree) -> bpy.types.GeometryNode:
    """
    Add a join last node to the tree.
    """
    link = tree.links.new
    node_join: bpy.types.GeometryNode = tree.nodes.new("GeometryNodeJoinGeometry")
    node_output = get_output(tree)
    old_loc = node_output.location.copy()
    node_output.location += Vector([NODE_SPACING * 2, 0])
    node_join.location = old_loc + Vector([NODE_SPACING, 0])
    try:
        if len(node_output.inputs[0].links) > 0:
            from_socket = node_output.inputs[0].links[0].from_socket
            if from_socket.node != get_input(tree):
                link(
                    node_output.inputs[0].links[0].from_socket,
                    node_join.inputs[0],
                )
    except IndexError:
        link(node_join.outputs[0], node_output.inputs[0])

    link(node_join.outputs[0], tree.nodes["Group Output"].inputs[0])
    return node_join


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
        # insert_join_last(tree)


def loc_between(a: bpy.types.GeometryNode, b: bpy.types.GeometryNode, t=0.5) -> Vector:
    """
    Get the location between two nodes
    """
    return a.location + (b.location - a.location) * t


def insert_before(
    item: bpy.types.Node | bpy.types.NodeSocket,
    new_node: str,
    offset: Vector = Vector([-NODE_SPACING, 0]),
) -> bpy.types.Node | bpy.types.GeometryNodeGroup:
    """
    Place a node before the given node in the tree. If a socket is given, link to that
    socket otherwise e link to the first input of the node.
    """
    # if a socket is given, then we will link into that socket, but if a node is given
    # we move down through the inputs and find the first one that is linked and link into
    # that socket
    if isinstance(item, bpy.types.NodeSocket):
        to_socket = item
        node = to_socket.node
        try:
            from_socket = to_socket.links[0].from_socket  # type: ignore
        except IndexError:
            from_socket = None
    else:
        node = item
        to_socket = node.inputs[0]
        from_socket = to_socket.links[0].from_socket

    tree = node.id_data
    try:
        node_new = add_custom(tree, new_node)
    except KeyError:
        node_new = tree.nodes.new(new_node)

    node_new.location = node.location + offset
    tree.links.new(node_new.outputs[0], to_socket)
    if from_socket is not None:
        tree.links.new(from_socket, node_new.inputs[0])
    return node_new


def custom_boolean_iswitch(
    name: str,
    items: Iterable[str],
    attribute_name: str = "chain_id",
    offset: int = 0,
    prefix: str = "",
) -> bpy.types.GeometryNodeTree:

    with g.tree(name) as tree:
        attr = g.NamedAttribute.integer(attribute_name)

        switch = g.IndexSwitch.boolean(
            index=attr if offset == 0 else attr + offset,
            items=[tree.inputs.boolean(prefix + x) for x in items],
        )

        switch >> tree.outputs.boolean("Selection")
        ~switch >> tree.outputs.boolean("Inverted")

    tree.tree.color_tag = "INPUT"

    return tree.tree


def custom_color_iswitch(
    name: str,
    items: dict[str, tuple[float, float, float, float]] | Iterable[int | float],
    attribute_name: str = "chain_id",
    offset: int = 0,
) -> bpy.types.GeometryNodeTree:
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
    return tree.tree
