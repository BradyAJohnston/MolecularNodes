import bpy
from bpy.types import Node
from mathutils import Vector
from typing import List, Sequence
import re

from . import nodes
from .utils import socket, option, TreeInterface
import numpy as np

NODE_SPACING = 250


def insert_join_last(tree: bpy.types.GeometryNodeTree) -> bpy.types.Node:
    """
    Add a join last node to the tree.
    """
    link = tree.links.new
    node_join = tree.nodes.new("GeometryNodeJoinGeometry")
    node_output = nodes.get_output(tree)
    old_loc = node_output.location.copy()
    node_output.location += Vector([NODE_SPACING * 2, 0])
    node_join.location = old_loc + Vector([NODE_SPACING, 0])
    try:
        if len(node_output.inputs[0].links) > 0:
            from_socket = node_output.inputs[0].links[0].from_socket  # type: ignore
            if from_socket.node != nodes.get_input(tree):
                link(
                    node_output.inputs[0].links[0].from_socket,  # type: ignore
                    node_join.inputs[0],
                )
    except IndexError:
        pass

    link(node_join.outputs[0], tree.nodes["Group Output"].inputs[0])
    return node_join


def final_join(tree: bpy.types.GeometryNodeTree) -> bpy.types.Node:
    """
    Get the last JoinGeometry node in the tree.
    """
    output = nodes.get_output(tree)
    try:
        linked = output.inputs[0].links[0].from_socket.node  # type: ignore
        if linked.bl_idname == "GeometryNodeJoinGeometry":
            return linked
        else:
            return insert_join_last(tree)
    except IndexError:
        return insert_join_last(tree)


def loc_between(a: bpy.types.Node, b: bpy.types.Node, t=0.5) -> Vector:
    """
    Get the location between two nodes
    """
    return a.location + (b.location - a.location) * t


def input_named_attribute(
    socket: bpy.types.NodeSocket, name: str, data_type: str | None
) -> bpy.types.Node:
    """
    Add a named attribute node to the tree and connect it to the given socket
    """
    tree = socket.node.id_data
    node_na = tree.nodes.new("GeometryNodeInputNamedAttribute")

    if data_type is not None:
        node_na.data_type = data_type
    node_na.inputs["Name"].default_value = name
    node_na.location = socket.node.location - Vector([NODE_SPACING, 0])

    tree.links.new(node_na.outputs["Attribute"], socket)

    return node_na


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
        node_new = nodes.add_custom(tree, new_node)
    except Exception:
        node_new = tree.nodes.new(node.bl_idname)

    node_new.location = node.location + offset
    tree.links.new(node_new.outputs[0], to_socket)
    if from_socket is not None:
        tree.links.new(from_socket, node_new.inputs[0])
    return node_new


def insert_set_color(
    node: bpy.types.Node,
    color: str | Sequence[float] = [0.3, 0.3, 0.3, 1.0],
) -> bpy.types.GeometryNodeGroup:
    """
    Add a set color node to the tree and connect it to the given socket
    """
    tree = node.id_data
    node_sc: bpy.types.GeometryNodeGroup = insert_before(node, "Set Color")  # type: ignore

    if color.lower() in ["default", "common"]:
        node_cc = insert_before(node_sc.inputs["Color"], "Color Common")
        node_car: bpy.types.GeometryNodeGroup = insert_before(  # type: ignore
            node_cc.inputs["Carbon"], "Color Attribute Random"
        )

        return node_car

    if isinstance(color, str):
        input_named_attribute(node_sc.inputs["Color"], color, "FLOAT_COLOR")
    else:
        node_sc.inputs["Color"].default_value = color
    return node_sc


def get_links(socket: bpy.types.NodeSocket) -> bpy.types.NodeLinks:
    """
    Get the links between two sockets
    """
    links = socket.links
    if links is None:
        raise ValueError("No links found for socket")
    return links


def add_style_branch(
    tree: bpy.types.GeometryNodeTree,
    style: str,
    color: str | None = None,
    selection: str | None = None,
) -> None:
    """
    Add a style branch to the tree.
    """
    style_name = nodes.styles_mapping[style]
    link = tree.links.new
    input = nodes.get_input(tree)
    output = nodes.get_output(tree)
    node_join = final_join(tree)

    current_min_y = min(node.location[1] for node in tree.nodes)
    ypos = current_min_y - 200
    xpos = loc_between(input, node_join, 0.75)[0]

    # ensure node exists
    nodes.append(style_name)
    node_style = nodes.add_custom(
        group=tree,
        name=style_name,
        location=[xpos, ypos],
    )

    link(
        input.outputs[0],
        node_style.inputs[0],
    )
    link(
        node_style.outputs[0],
        node_join.inputs[0],
    )
    if selection:
        input_named_attribute(node_style.inputs["Selection"], selection, "BOOLEAN")
    if color:
        insert_set_color(node_style, color)


def get_final_style_nodes(tree: bpy.types.GeometryNodeTree) -> List[bpy.types.Node]:
    """
    Get the final style nodes in the tree.
    """
    links: bpy.types.NodeLinks = final_join(tree).inputs[0].links  # type: ignore
    return [
        link.from_socket.node
        for link in links
        if link.from_socket.node.name.startswith("Style")
    ]


class GeometryNodeInterFace(TreeInterface):
    """
    Interface for the geometry nodes in the tree.
    """

    def __init__(self, node: bpy.types.Node) -> None:
        self.tree = node.id_data

    @property
    def node_tree(self) -> bpy.types.GeometryNodeTree:
        return self.tree

    def _expose_options(self, node: bpy.types.Node) -> None:
        for input in node.inputs:
            if input.is_linked:
                continue
            try:
                value = input.default_value  # type: ignore
            except AttributeError:
                continue

            # Determine the type based on the value
            if isinstance(value, (int, float, bool, str)):
                value_type = type(value)
            else:
                value_type = np.ndarray
                value = np.array(value)

            # Create property name from node and input names
            prop_name = (
                "_".join([node.name.split(".")[0], input.name])
                .lower()
                .replace(" ", "_")
                .replace("style_", "")
            )
            # Add the socket as a property to the class
            # shader sockets aren't to be accessed directly so just ignore them
            try:
                prop = socket(node.name, input.name, value_type)  # type: ignore
                setattr(self.__class__, prop_name, prop)
            except AttributeError:
                pass


def create_style_interface(node: Node, linked: bool = True) -> GeometryNodeInterFace:
    """
    Dynamically create a StyleInterface class with exposed options.
    """
    class_name = f"DynamicStyleInterface_{node.name}"

    # Create the dynamic class with a unique name
    DynamicInterface = type(class_name, (GeometryNodeInterFace,), {})

    # Pre-expose the options for the main node
    interface = DynamicInterface(node)
    interface._expose_options(node)

    if not linked:
        return interface

    # Pre-expose options for linked nodes
    for input in node.inputs:
        if input.is_linked:
            linked_node = input.links[0].from_socket.node  # type: ignore
            interface._expose_options(linked_node)

    return interface


class StyleWrangler:
    """
    Class to manage the style nodes in the tree.
    """

    def __init__(self, tree: bpy.types.GeometryNodeTree):
        self.tree = tree

    @property
    def styles(self) -> List[GeometryNodeInterFace]:
        """
        Get the styles in the tree.
        """
        return [
            create_style_interface(node) for node in get_final_style_nodes(self.tree)
        ]
