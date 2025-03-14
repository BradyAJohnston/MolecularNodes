import bpy
from mathutils import Vector
from typing import List

from . import nodes


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
    color: str | List[float, float, float, float],
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
