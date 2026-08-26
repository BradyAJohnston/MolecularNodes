from typing import List
import bpy
from bpy.types import Node
from databpy.nodes import get_output
from .arrange import arrange_tree
from .nodes import (
    final_join,
)


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
