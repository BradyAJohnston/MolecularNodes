from typing import List, Sequence
import bpy
from bpy.types import Node
from databpy.nodes import get_output
from mathutils import Vector
from . import nodes
from .arrange import arrange_tree
from .interface import (
    input_named_attribute,
)
from .nodes import (
    NODE_SPACING,
    final_join,
    insert_before,
)


def insert_set_color(
    node: bpy.types.Node,
    color: str | Sequence[float] = [0.3, 0.3, 0.3, 1.0],
) -> bpy.types.GeometryNodeGroup:
    """
    Add a set color node to the tree and connect it to the given socket
    """
    _tree = node.id_data
    node_sc: bpy.types.GeometryNodeGroup = insert_before(node, "Set Color")  # type: ignore

    if isinstance(color, str):
        if color.lower() in ["default", "common"]:
            node_cc = insert_before(node_sc.inputs["Color"], "Color Common")
            node_car: bpy.types.GeometryNodeGroup = insert_before(  # type: ignore
                node_cc.inputs["Carbon"], "Color Attribute Random"
            )
            return node_car
        elif color.lower() == "plddt":
            insert_before(node_sc.inputs["Color"], "Color pLDDT")
        else:
            input_named_attribute(node_sc.inputs["Color"], color, "FLOAT_COLOR")
    else:
        node_sc.inputs["Color"].default_value = color
    return node_sc


def insert_animate_frames(
    node: bpy.types.GeometryNode, frames: bpy.types.Collection | str
) -> bpy.types.Node:
    """
    Add an animate frames node to the tree and connect it to the given socket
    """
    node.location += Vector([NODE_SPACING, 0])
    tree = node.id_data
    node_af: bpy.types.GeometryNodeGroup = insert_before(node, "Animate Frames")  # type: ignore
    if isinstance(frames, bpy.types.Collection):
        node_af.inputs["Frames"].default_value = frames
    elif isinstance(frames, str):
        node_af.inputs["Frames"].default_value = bpy.data.collections[frames]
    else:
        raise ValueError(
            f"Frames must be a string or a Collection, not {type(frames)=}"
        )

    node_an = nodes.add_custom(
        tree,
        "Animate Value",
        location=node_af.location + Vector([-NODE_SPACING, -NODE_SPACING / 1.5]),
    )
    node_an.inputs["Value Max"].default_value = float(len(frames.objects) - 1)  # type: ignore
    tree.links.new(node_an.outputs[0], node_af.inputs["Frame"])

    return node_af


def get_final_style_nodes(
    tree: bpy.types.GeometryNodeTree,
) -> List[bpy.types.GeometryNodeGroup]:
    """
    Get the final style nodes in the tree.
    """
    try:
        links: bpy.types.NodeLinks = final_join(tree).inputs[0].links
    except AttributeError:
        links = get_output(tree).inputs[0].links

    return [
        link.from_socket.node
        for link in reversed(links)
        if link.from_socket.node.name.startswith("Style")
    ]


def remove_style_node(node: Node) -> None:
    """
    Remove a style node from its tree, along with the nodes linked into its inputs.
    """
    tree: bpy.types.NodeTree = node.id_data
    to_remove = [node] + [
        input.links[0].from_socket.node for input in node.inputs if input.is_linked
    ]
    for node_to_remove in to_remove:
        tree.nodes.remove(node_to_remove)
    arrange_tree(tree)
