import bpy
from bpy.types import Node
from mathutils import Vector
from typing import List
import re

from . import nodes
from .utils import socket, option, TreeInterface
import numpy as np


def insert_join_last(tree: bpy.types.GeometryNodeTree) -> bpy.types.Node:
    """
    Add a join last node to the tree.
    """
    node_join = tree.nodes.new("GeometryNodeJoinGeometry")
    node_join.location = [
        tree.nodes["Group Output"].location[0] - 200,
        tree.nodes["Group Output"].location[1],
    ]
    link = tree.links.new
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


def add_style_branch(
    tree: bpy.types.GeometryNodeTree,
    style: str,
    color: str | None = None,
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

    def __init__(self, tree: bpy.types.GeometryNodeTree) -> None:
        self.tree = tree

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
                setattr(
                    self,
                    prop_name,
                    socket(node.name, input.name, value_type),  # type: ignore
                )
            except AttributeError:
                pass


class StyleInterface(GeometryNodeInterFace):
    """
    Interface for the style nodes in the tree.
    """

    def __init__(self, node: Node):
        super().__init__(node.id_data)
        self.node = node

        # for an intial test we expose the options of the given node and also any nodes
        # that are directly linked into it. Could walk back more, but for testing now
        # we stick with just 1 level
        self._expose_options(node)
        for input in self.node.inputs:
            if input.is_linked:
                node = input.links[0].from_socket.node  # type: ignore
                self._expose_options(node)

    def __repr__(self) -> str:
        return f"<Interface for {self.node.name} in {self.node.id_data.name}>"


class StyleWrangler:
    """
    Class to manage the style nodes in the tree.
    """

    def __init__(self, tree: bpy.types.GeometryNodeTree):
        self.tree = tree

    @property
    def styles(self) -> List[StyleInterface]:
        """
        Get the styles in the tree.
        """
        return [StyleInterface(node) for node in get_final_style_nodes(self.tree)]
