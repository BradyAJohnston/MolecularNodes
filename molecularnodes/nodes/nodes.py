import math
from typing import Iterable
import bpy
import databpy
import numpy as np
from databpy.nodes import (
    append_from_blend,
)
from mathutils import Vector
from nodebpy import geometry as g
from .. import color, utils
from ..assets import MN_DATA_FILE
from ..blender import mesh
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


def outputs(node):
    items = {}
    for item in node.interface.items_tree:
        if item.item_type == "SOCKET":
            if item.in_out == "OUTPUT":
                items[item.name] = item
    return items


def get_output_type(node, type="INT"):
    for output in node.outputs:
        if output.type == type:
            return output


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


def get_color_node(object):
    "Walk back through the primary node connections until you find the first style node"
    group = object.modifiers["Molecular Nodes"].node_group
    for node in group.nodes:
        if node.name == "Color Attribute Random":
            return node


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
    )  # type: ignore
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
    node: bpy.types.GeometryNodeGroup = group.nodes.new("GeometryNodeGroup")  # type: ignore
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


def create_starting_node_tree(
    object: bpy.types.Object,
    coll_frames: bpy.types.Collection | None = None,
    style: str = "spheres",
    name: str | None = None,
    color: str | None = "common",
    material: str | bpy.types.Material = "MN Default",
    is_modifier: bool = True,
):
    """
    Create a starting node tree for the inputted object.

    Parameters
    ----------
    object : bpy.types.Object
        Object to create the node tree for.
    coll_frames : bpy.data.collections, optional
        If None, no animation will be created.
        The default is None.
    style : str, optional
        Starting style for the node tree. The default is "spheres".
        Available options are stored as the keys of styles_mapping.
    name : str, optional
        Name of the node tree. If None, a default name will be generated.
        The default is None.
    color : str, optional
        None doesn't add ay set_color nodes, 'common' adds the color by common elements
        and 'plddt' adds color by pLDDT score.
    """

    # ensure there is a geometry nodes modifier called 'MolecularNodes' that is created and applied to the object
    mod: bpy.types.NodesModifier = get_mod(object)

    if not name:
        name = f"MN_{object.name}"

    with g.tree(name) as tree:
        atoms = tree.inputs.geometry("Atoms")
        join = g.JoinGeometry()
        join >> tree.outputs.geometry("Geometry")

        match color.lower():
            case "pldtt":
                color = g.ColorPLDDT()
            case _:
                color = g.ColorElement(c=g.RandomColor(g.ChainID()))

        if coll_frames:
            atoms = atoms >> g.AnimateFrames(
                collection=coll_frames, factor=g.AnimateValue()
            )

        style_node = {
            "ribbon": g.StyleRibbon,
            "cartoon": g.StyleCartoon,
            "surface": g.StyleSurface,
            "ball_and_stick": g.StyleBallAndStick,
            "sticks": g.StyleSticks,
        }[style]

        if isinstance(material, str):
            material = bpy.data.materials()

        assign_material(style_node.node, material)

        atoms >> g.SetColor(color=color) >> style_node >> join

    mod.node_group = tree.tree


def combine_join_geometry(group, node_list, output="Geometry", join_offset=300):
    link = group.links.new
    max_x = max([node.location[0] for node in node_list])
    node_to_instances = group.nodes.new("GeometryNodeJoinGeometry")
    node_to_instances.location = [int(max_x + join_offset), 0]

    for node in reversed(node_list):
        link(node.outputs[output], node_to_instances.inputs["Geometry"])
    return node_to_instances


def split_geometry_to_instances(name, iter_list=("A", "B", "C"), attribute="chain_id"):
    """Create a Node to Split Geometry by an Attribute into Instances

    Splits the inputted geometry into instances, based on an attribute field. By
    default this field is the `chain_id` but this can be selected for any field.
    Will loop over each item of the list, so a list of arbitrary items that will
    define how many times to create the required nodes.

    """
    group = new_tree(name)
    node_input = get_input(group)
    node_output = get_output(group)

    named_att = group.nodes.new("GeometryNodeInputNamedAttribute")
    named_att.location = [-200, -200]
    named_att.data_type = "INT"
    named_att.inputs[0].default_value = attribute

    link = group.links.new
    list_sep = []

    for i, chain in enumerate(iter_list):
        pos = [i % 10, math.floor(i / 10)]

        node_split = add_custom(group, ".MN_utils_split_instance")
        node_split.location = [int(250 * pos[0]), int(-300 * pos[1])]
        node_split.inputs["Group ID"].default_value = i
        link(named_att.outputs["Attribute"], node_split.inputs["Field"])
        link(node_input.outputs[0], node_split.inputs["Geometry"])
        list_sep.append(node_split)

    node_instance = combine_join_geometry(group, list_sep, "Instance")
    node_output.location = [int(10 * 250 + 400), 0]
    link(node_instance.outputs[0], node_output.inputs[0])
    return group


def assembly_data_object_from_obj(obj: bpy.types.Object) -> bpy.types.Object:
    data_obj_name = f".data_{obj.name}_assemblies"
    data_obj = bpy.data.objects.get(data_obj_name)
    if not data_obj:
        transforms = utils.array_quaternions_from_dict(obj.mn.biological_assemblies)
        data_obj = mesh.create_data_object(array=transforms, name=data_obj_name)

    return data_obj


def assembly_initialise(obj: bpy.types.Object):
    """
    Setup the required data object and nodes for building an assembly.
    """
    data_obj = assembly_data_object_from_obj(obj)
    tree_assembly = create_assembly_node_tree(name=obj.name, data_object=data_obj)
    return tree_assembly


def assembly_insert(mol: bpy.types.Object):
    """
    Given a molecule, setup the required assembly node and insert it into the node tree.
    """

    tree_assembly = assembly_initialise(mol)
    group = get_mod(mol).node_group
    node = add_custom(group, tree_assembly.name)
    insert_last_node(get_mod(mol).node_group, node)


def create_assembly_node_tree(
    name: str, data_object: bpy.types.Object
) -> bpy.types.NodeTree:
    node_group_name = f"Assembly {name}"
    existing_node_tree = bpy.data.node_groups.get(node_group_name)
    if existing_node_tree:
        return existing_node_tree

    tree: bpy.types.NodeTree = new_tree(name=node_group_name)
    link = tree.links.new

    node_split = add_custom(tree, "Split to Centred Instances", [-150, 0])

    node_att: bpy.types.GeometryNodeInputNamedAttribute = tree.nodes.new(
        "GeometryNodeInputNamedAttribute"
    )  # type: ignore
    node_att.data_type = "INT"
    node_att.inputs[0].default_value = "chain_id"
    node_att.location = [-150, -200]
    link(node_att.outputs["Attribute"], node_split.inputs["Group ID"])

    node_group_assembly_instance = append(".MN_assembly_instance_chains")
    node_assembly = add_custom(tree, node_group_assembly_instance.name, [150, 0])
    node_assembly.inputs["data_object"].default_value = data_object

    out_sockets = outputs(tree)
    out_sockets[list(out_sockets)[0]].name = "Instances"

    socket_info = (
        {
            "name": "Rotation",
            "type": "NodeSocketFloat",
            "min": 0,
            "max": 1,
            "default": 1,
        },
        {
            "name": "Translation",
            "type": "NodeSocketFloat",
            "min": 0,
            "max": 1,
            "default": 1,
        },
        {
            "name": "assembly_id",
            "type": "NodeSocketInt",
            "min": 1,
            "max": max(databpy.named_attribute(data_object, "assembly_id")),
            "default": 1,
        },
    )

    for info in socket_info:
        socket = tree.interface.items_tree.get(info["name"])
        if not socket:
            socket: bpy.types.NodeTreeInterfaceSocket = tree.interface.new_socket(
                info["name"], in_out="INPUT", socket_type=info["type"]
            )
        socket.default_value = info["default"]
        socket.min_value = info["min"]
        socket.max_value = info["max"]

        link(get_input(tree).outputs[info["name"]], node_assembly.inputs[info["name"]])

    get_output(tree).location = [400, 0]
    link(get_input(tree).outputs[0], node_split.inputs[0])
    link(node_split.outputs[0], node_assembly.inputs[0])
    link(node_assembly.outputs[0], get_output(tree).inputs[0])
    if hasattr(tree, "color_tag"):
        tree.color_tag = "GEOMETRY"
    return tree


def add_inverse_selection(group):
    output = get_output(group)
    if "Inverted" not in output.inputs.keys():
        group.interface.new_socket(
            "Inverted", in_out="OUTPUT", socket_type="NodeSocketBool"
        )

    loc = output.location
    bool_math = group.nodes.new("FunctionNodeBooleanMath")
    bool_math.location = [loc[0], -100]
    bool_math.operation = "NOT"

    group.links.new(
        output.inputs["Selection"].links[0].from_socket, bool_math.inputs[0]
    )
    group.links.new(bool_math.outputs[0], output.inputs["Inverted"])


def boolean_link_output(tree: bpy.types.NodeTree, node: bpy.types.Node) -> None:
    link = tree.links.new
    node_output = get_output(tree)
    tree.interface.new_socket(
        name="Selection", in_out="OUTPUT", socket_type=socket_types["BOOLEAN"]
    )
    tree.interface.new_socket(
        name="Inverted", in_out="OUTPUT", socket_type=socket_types["BOOLEAN"]
    )
    final_output = node.outputs[0]
    link(final_output, node_output.inputs["Selection"])
    node_invert: bpy.types.FunctionNodeBooleanMath = tree.nodes.new(
        "FunctionNodeBooleanMath"
    )  # type: ignore

    node_invert.operation = "NOT"
    node_invert.location = (np.array(node_output.location) - [0, 200]).tolist()
    link(final_output, node_invert.inputs[0])
    link(node_invert.outputs[0], node_output.inputs["Inverted"])


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
        # for socket in node.inputs:
        #     if socket.is_linked:
        #         from_socket = socket.links[0].from_socket
        #         to_socket = socket
        #         break

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
