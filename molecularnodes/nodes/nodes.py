import itertools
import math
from typing import List, Optional
import bpy
import databpy
import numpy as np
from databpy.nodes import (
    NodeGroupCreationError,
    append_from_blend,
    swap_tree,
)
from mathutils import Vector
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
    "density_surface": "Style Density Surface",
    "density_wire": "Style Density Wire",
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
        custom_iswitch(
            name="selection", iter_list=input_list, field=field, dtype="BOOLEAN"
        ).name,
    )

    set_selection(group, style, sel_node)
    return sel_node


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


def get_mod(object, name="MolecularNodes"):
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


def style_node(group):
    prev = previous_node(get_output(group))
    is_style_node = "Style" in prev.name
    while not is_style_node:
        prev = previous_node(prev)
        is_style_node = "Style" in prev.name
    return prev


def get_style_node(object):
    "Walk back through the primary node connections until you find the first style node"
    group = object.modifiers["MolecularNodes"].node_group
    return style_node(group)


def star_node(group):
    prev = previous_node(get_output(group))
    is_star_node = "Starfile Instances" in prev.name
    while not is_star_node:
        prev = previous_node(prev)
        is_star_node = "Starfile Instances" in prev.name
    return prev


def get_star_node(object):
    "Walk back through the primary node connections until you find the first style node"
    group = object.modifiers["MolecularNodes"].node_group
    return star_node(group)


def get_color_node(object):
    "Walk back through the primary node connections until you find the first style node"
    group = object.modifiers["MolecularNodes"].node_group
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
    group = obj.modifiers["MolecularNodes"].node_group
    realize = group.nodes.new("GeometryNodeRealizeInstances")
    insert_last_node(group, realize)


def swap(node: bpy.types.Node, tree: str | bpy.types.NodeTree) -> None:
    "Swap out the node's node_tree, while maintaining the possible old connections"

    if isinstance(tree, str):
        try:
            tree = bpy.data.node_groups[tree]
        except KeyError:
            tree = append(tree)

    swap_tree(node=node, tree=tree)


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
    tree = bpy.data.node_groups.get(name)  # type: ignore
    # if the group already exists, return it and don't create a new one
    if tree and fallback:
        if not isinstance(tree, bpy.types.GeometryNodeTree):
            raise TypeError(f"Expected a GeometryNodeTree, got {type(tree)}")
        return tree

    # create a new group for this particular name and do some initial setup
    tree: bpy.types.GeometryNodeTree = bpy.data.node_groups.new(
        name=name,
        type="GeometryNodeTree",  # type: ignore
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

    # if there is an input socket called 'Material', assign it to the base MN material
    # if another material is supplied, use that instead.
    assign_material(node, new_material=material)

    # move and format the node for arranging
    node.location = location
    node.width = width
    node.show_options = show_options
    node.name = name

    return node


def change_style_node(obj: bpy.types.Object, style: str):
    swap(get_style_node(obj), append(styles_mapping[style]))


def create_starting_nodes_starfile(object):
    # ensure there is a geometry nodes modifier called 'MolecularNodes' that is created and applied to the object
    node_mod = get_mod(object)

    node_name = f"MN_starfile_{object.name}"

    # create a new GN node group, specific to this particular molecule
    group = new_tree(node_name)
    node_mod.node_group = group
    link = group.links.new

    # move the input and output nodes for the group
    node_input = get_input(group)
    node_output = get_output(group)
    node_input.location = [0, 0]
    node_output.location = [700, 0]
    node_star_instances = add_custom(group, "Starfile Instances", [450, 0])
    link(node_star_instances.outputs[0], node_output.inputs[0])
    link(node_input.outputs[0], node_star_instances.inputs[0])


def create_starting_nodes_density(object, threshold=0.8, style="density_surface"):
    # ensure there is a geometry nodes modifier called 'MolecularNodes' that is created and applied to the object
    mod = get_mod(object)
    node_name = f"MN_density_{object.name}"

    try:
        tree = bpy.data.node_groups[node_name]
        mod.node_group = tree
        return
    except KeyError:
        pass

    # create a new GN node group, specific to this particular molecule
    group = new_tree(node_name, fallback=False)
    link = group.links.new
    mod.node_group = group

    # move the input and output nodes for the group
    node_input = get_input(group)
    node_input.location = [0, 0]
    node_output = get_output(group)
    node_output.location = [800, 0]

    node_density = add_custom(group, styles_mapping[style], [400, 0])
    node_density.inputs["Threshold"].default_value = threshold

    link(node_input.outputs[0], node_density.inputs[0])
    link(node_density.outputs[0], node_output.inputs[0])


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
    mod = get_mod(object)

    if not name:
        name = f"MN_{object.name}"

    # check if the node tree already exists and use that instead
    # try:
    #     tree = bpy.data.node_groups[name]
    #     mod.node_group = tree
    #     return
    # except KeyError:
    #     pass

    tree = new_tree(name, input_name="Atoms")
    tree.is_modifier = is_modifier
    link = tree.links.new
    mod.node_group = tree

    # move the input and output nodes for the group
    node_input = get_input(tree)
    node_output = get_output(tree)
    node_input.location = [0, 0]
    node_output.location = [700, 0]

    if style is None:
        link(node_input.outputs[0], node_output.inputs[0])
        return tree

    node_style = add_custom(tree, styles_mapping[style], [450, 0], material=material)
    link(node_style.outputs[0], node_output.inputs[0])
    link(node_input.outputs[0], node_style.inputs[0])

    # if requested, setup the nodes for generating colors in the node tree
    if color is not None:
        if color == "common":
            node_color_set = add_custom(tree, "Set Color", [200, 0])
            node_color_common = add_custom(tree, "Color Common", [-50, -150])
            node_random_color = add_custom(tree, "Color Attribute Random", [-300, -150])

            link(node_input.outputs[0], node_color_set.inputs[0])
            link(node_random_color.outputs["Color"], node_color_common.inputs["Carbon"])
            link(node_color_common.outputs[0], node_color_set.inputs["Color"])
            link(node_color_set.outputs[0], node_style.inputs[0])
            to_animate = node_color_set
        elif color.lower() == "plddt":
            node_color_set = add_custom(tree, "Set Color", [200, 0])
            node_color_plddt = add_custom(tree, "Color pLDDT", [-50, -150])

            link(node_input.outputs[0], node_color_set.inputs["Atoms"])
            link(node_color_plddt.outputs[0], node_color_set.inputs["Color"])
            link(node_color_set.outputs["Atoms"], node_style.inputs["Atoms"])
        else:
            to_animate = node_input
    else:
        to_animate = node_input

    # if multiple frames, set up the required nodes for an animation
    if coll_frames:
        node_output.location = [1100, 0]
        node_style.location = [800, 0]

        node_animate_frames = add_custom(tree, "Animate Frames", [500, 0])
        node_animate = add_custom(tree, "Animate Value", [500, -300])

        node_animate_frames.inputs["Frames"].default_value = coll_frames
        node_animate.inputs["Value Max"].default_value = len(coll_frames.objects) - 1

        link(to_animate.outputs[0], node_animate_frames.inputs[0])
        link(node_animate_frames.outputs[0], node_style.inputs[0])
        link(node_animate.outputs[0], node_animate_frames.inputs["Frame"])


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


def assembly_initialise(obj: bpy.types.Object):
    """
    Setup the required data object and nodes for building an assembly.
    """

    data_obj_name = f".data_assembly_{obj.name}"

    # check if a data object exists and create a new one if not
    data_object = bpy.data.objects.get(data_obj_name)
    if not data_object:
        transforms = utils.array_quaternions_from_dict(obj.mn.biological_assemblies)
        data_object = mesh.create_data_object(array=transforms, name=data_obj_name)

    tree_assembly = create_assembly_node_tree(name=obj.name, data_object=data_object)
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
    node_join: bpy.types.GeometryNode = tree.nodes.new("GeometryNodeJoinGeometry")  # type: ignore
    node_output = get_output(tree)
    old_loc = node_output.location.copy()
    node_output.location += Vector([NODE_SPACING * 2, 0])
    node_join.location = old_loc + Vector([NODE_SPACING, 0])
    try:
        if len(node_output.inputs[0].links) > 0:
            from_socket = node_output.inputs[0].links[0].from_socket  # type: ignore
            if from_socket.node != get_input(tree):
                link(
                    node_output.inputs[0].links[0].from_socket,  # type: ignore
                    node_join.inputs[0],
                )
    except IndexError:
        pass

    link(node_join.outputs[0], tree.nodes["Group Output"].inputs[0])
    return node_join


def final_join(tree: bpy.types.GeometryNodeTree) -> bpy.types.GeometryNode:
    """
    Get the last JoinGeometry node in the tree.
    """
    output = get_output(tree)
    if "Assembly" in output.inputs[0].links[0].from_socket.node.name:
        if (
            output.inputs[0].from_socket.node.inputs[0].links[0].from_socket.node.name
            == "GeometryNodeJoinGeometry"
        ):
            return output.inputs[0].from_socket.node.inputs[0].from_socket.node  # type: ignore
    try:
        linked = output.inputs[0].links[0].from_socket.node  # type: ignore
        if linked.bl_idname == "GeometryNodeJoinGeometry":
            return linked
        else:
            return insert_join_last(tree)
    except IndexError:
        return insert_join_last(tree)


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
        from_socket = to_socket.links[0].from_socket  # type: ignore
        # for socket in node.inputs:
        #     if socket.is_linked:
        #         from_socket = socket.links[0].from_socket  # type: ignore
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


def custom_iswitch(
    name: str,
    iter_list,
    field: str = "chain_id",
    dtype: str = "BOOLEAN",
    default_values=None,
    prefix: str = "",
    start: int = 0,
    offset: int = 0,
    panels: Optional[List[str]] = None,
    panels_open: int = 1,
) -> bpy.types.GeometryNodeTree:
    """
    Creates a named `Index Switch` node.

    Wraps an index switch node, giving the group names or each name in the `iter_list`. The
    inputs can also be placed in subpanels and given specific default values. Uses the
    given field for the attribute name to use in the index switch, and optionally adds an
    offset value if the start value is non zero.

    If a list of default items is given, then it is recycled to fill the defaults for
    each created socket in for the node.

    Parameters
    ----------
    name : str
        The name of the node group.
    iter_list : list
        The list of items to iterate over.
    field : str, optional
        The name of the attribute field. Defaults to 'chain_id'.
    default_values : list, optional
        The list of default values to assign to each item. Defaults to None.
    panels : list, str
        List of panel names for the sockets to be assigned to. If None, then socket will
        not be assigned to the panel. The will appear in the panel in the order in which
        they are given.
    panels_open: int
        Number of panels to default to open. If `0` then all panels will be closed. Useful
        for larger panel node groups, to keep them closed and help with organisation. A
        list of panel names can also be given to optionally categorise the sockets into
        panels in the final group node. The length of the panels list must match the length
        of the iter_list. `None` values mean the socket is not placed in a panel.
    prefix : str, optional
        The prefix to add to the node names. Defaults to an empty string.
    start : int, optional
        The starting index for the node names. Defaults to 0.

    Returns
    -------
    group : bpy.types.GeometryNodeTree
        The created node group.

    Raises
    ------
    NodeGroupCreationError
        If there was an error creating the node group.
    """
    iter_list = [str(i) for i in iter_list]
    tree = bpy.data.node_groups.get(name)
    if tree:
        return tree

    socket_type = socket_types[dtype]
    tree = new_tree(name, geometry=False, fallback=False)

    # try creating the node group, otherwise on fail cleanup the created group and
    # report the error
    try:
        link = tree.links.new
        node_input = get_input(tree)
        node_output = get_output(tree)
        node_attr = tree.nodes.new("GeometryNodeInputNamedAttribute")
        node_attr.data_type = "INT"
        node_attr.location = [0, 150]
        node_attr.inputs["Name"].default_value = str(field)

        node_iswitch: bpy.types.GeometryNodeIndexSwitch = tree.nodes.new(  # type: ignore
            "GeometryNodeIndexSwitch"
        )
        node_iswitch.data_type = dtype

        link(node_attr.outputs["Attribute"], node_iswitch.inputs["Index"])

        # if there is as offset to the lookup values (say we want to start looking up
        # from 100 or 1000 etc) then we add a math node with that offset value
        if start != 0:
            node_math = tree.nodes.new("ShaderNodeMath")
            node_math.operation = "ADD"
            node_math.location = [0, 150]
            node_attr.location = [0, 300]

            node_math.inputs[1].default_value = start
            link(node_attr.outputs["Attribute"], node_math.inputs[0])
            link(node_math.outputs["Value"], node_iswitch.inputs["Index"])

        # if there are custom values provided, create a dictionary lookup for those values
        # to assign to the sockets upon creation. If no default was given and the dtype
        # is colors, then generate a random pastel color for each value
        default_lookup = None
        if default_values is not None:
            default_lookup = dict(zip(iter_list, itertools.cycle(default_values)))
        elif dtype == "RGBA":
            default_lookup = dict(
                zip(iter_list, [color.random_rgb() for i in iter_list])
            )

        # for each item in the iter_list, we create a new socket on the interface for this
        # node group, and link it to the interface on the index switch. The index switch
        # currently starts with two items already, so once i > 1 we start to add
        # new items for the index switch as well
        panel_item_counter = 0
        panel_counter = 0
        for j in range(offset):
            node_iswitch.index_switch_items.new()
        for i, item in enumerate(iter_list):
            if i > 1:
                node_iswitch.index_switch_items.new()

            # The offset creates but skips itmes on the index switch node.
            # if we offset by 1 then we can index from 1 essentially, for attributes like
            # the atomic_number

            socket = tree.interface.new_socket(
                name=f"{prefix}{item}", in_out="INPUT", socket_type=socket_type
            )
            #  if a set of default values was given, then use it for setting
            # the defaults on the created sockets of the node group
            if default_lookup:
                socket.default_value = default_lookup[item]
            link(
                node_input.outputs[socket.identifier],
                node_iswitch.inputs[str(i + offset)],
            )

            # if a list of panel names has been passed in, then we use it to
            # assign all of the interface sockets to the panels
            if panels:
                pname = panels[i]
                if pname is not None:
                    # try and get an existing panel, if None is returned we have to create
                    # a panel with the given name
                    panel: bpy.types.NodeTreeInterfacePanel = (
                        tree.interface.items_tree.get(pname)
                    )
                    if not panel:
                        panel_item_counter = 0
                        panel = tree.interface.new_panel(name=pname)
                        panel_counter += 1
                        # we can set a certain number of panels to be open when created.
                        # a value of 0 means all created panels will be closed on creation.
                        # larger cutoffs will mean that n number of panels will be open by
                        # default
                        if panel_counter >= panels_open:
                            panel.default_closed = True

                    tree.interface.move_to_parent(
                        socket, panel, to_position=panel_item_counter + 1
                    )
                    panel_item_counter += 1

        if dtype == "BOOLEAN":
            tree.color_tag = "INPUT"
            boolean_link_output(tree, node_iswitch)
        elif dtype == "RGBA":
            tree.color_tag = "COLOR"
            socket_out = tree.interface.new_socket(
                name="Color", in_out="OUTPUT", socket_type=socket_type
            )
            link(
                node_iswitch.outputs["Output"],
                node_output.inputs[socket_out.identifier],
            )
        else:
            raise ValueError(f"Unsupported value typee for custom iswitch: {dtype}")

        return tree

    # if something broke when creating the node group, delete whatever was created
    except Exception as e:
        node_name = tree.name
        bpy.data.node_groups.remove(tree)
        raise NodeGroupCreationError(
            f"Unable to make node group: {node_name}.\nError: {e}"
        )


def resid_multiple_selection(node_name, input_resid_string):
    """
    Returns a node group that takes an integer input and creates a boolean
    tick box for each item in the input list. Outputs are the selected
    residues and the inverse selection. Used for constructing chain
    selections in specific proteins.
    """

    # do a cleanning of input string to allow fuzzy input from users
    for c in ";/+ .":
        if c in input_resid_string:
            input_resid_string = input_resid_string.replace(c, ",")

    for c in "_=:":
        if c in input_resid_string:
            input_resid_string = input_resid_string.replace(c, "-")

    # parse input_resid_string into sub selecting string list
    sub_list = [item for item in input_resid_string.split(",") if item]

    # distance vertical to space all of the created nodes
    node_sep_dis = -100

    # get the active object, might need to change to taking an object as an input
    # and making it active isntead, to be more readily applied to multiple objects

    # create the custom node group data block, where everything will go
    # also create the required group node input and position it
    residue_id_group = bpy.data.node_groups.new(node_name, "GeometryNodeTree")
    node_input = residue_id_group.nodes.new("NodeGroupInput")
    node_input.location = [0, node_sep_dis * len(sub_list) / 2]

    group_link = residue_id_group.links.new
    new_node = residue_id_group.nodes.new

    prev = None
    for residue_id_index, residue_id in enumerate(sub_list):
        # add an new node of Select Res ID or MN_sek_res_id_range
        current_node = new_node("GeometryNodeGroup")

        # add an bool_math block
        bool_math = new_node("FunctionNodeBooleanMath")
        bool_math.location = [400, (residue_id_index + 1) * node_sep_dis]
        bool_math.operation = "OR"

        if "-" in residue_id:
            # selecting a range of residues by using the Res ID Range node and connecting
            # to the min and max of those nodes
            current_node.node_tree = append("Select Res ID Range")
            [resid_start, resid_end] = residue_id.split("-")[:2]
            socket_1 = residue_id_group.interface.new_socket(
                "res_id: Min", in_out="INPUT", socket_type="NodeSocketInt"
            )
            socket_1.default_value = int(resid_start)
            socket_2 = residue_id_group.interface.new_socket(
                "res_id: Max", in_out="INPUT", socket_type="NodeSocketInt"
            )
            socket_2.default_value = int(resid_end)

            # a residue range
            group_link(
                node_input.outputs[socket_1.identifier], current_node.inputs["Min"]
            )
            group_link(
                node_input.outputs[socket_2.identifier], current_node.inputs["Max"]
            )
        else:
            # Selecting singular res ID numbers by creating the socket and adding a node
            # ensuring that we are connecting to the right node
            current_node.node_tree = append("Select Res ID")
            socket = residue_id_group.interface.new_socket(
                "res_id", in_out="INPUT", socket_type="NodeSocketInt"
            )
            socket.default_value = int(residue_id)
            group_link(
                node_input.outputs[socket.identifier], current_node.inputs["Res ID"]
            )

        # set the coordinates
        current_node.location = [200, (residue_id_index + 1) * node_sep_dis]
        if not prev:
            # link the first residue selection to the first input of its OR block
            group_link(current_node.outputs["Selection"], bool_math.inputs[0])
        else:
            # if it is not the first residue selection, link the output to the previous or block
            group_link(current_node.outputs["Selection"], prev.inputs[1])

            # link the ouput of previous OR block to the current OR block
            group_link(prev.outputs[0], bool_math.inputs[0])
        prev = bool_math

    # add a output block
    residue_id_group_out = new_node("NodeGroupOutput")
    residue_id_group_out.location = [800, (residue_id_index + 1) / 2 * node_sep_dis]
    residue_id_group.interface.new_socket(
        "Selection", in_out="OUTPUT", socket_type="NodeSocketBool"
    )
    residue_id_group.interface.new_socket(
        "Inverted", in_out="OUTPUT", socket_type="NodeSocketBool"
    )
    group_link(prev.outputs[0], residue_id_group_out.inputs["Selection"])
    invert_bool_math = new_node("FunctionNodeBooleanMath")
    invert_bool_math.location = [600, (residue_id_index + 1) / 3 * 2 * node_sep_dis]
    invert_bool_math.operation = "NOT"
    group_link(prev.outputs[0], invert_bool_math.inputs[0])
    group_link(invert_bool_math.outputs[0], residue_id_group_out.inputs["Inverted"])
    return residue_id_group
