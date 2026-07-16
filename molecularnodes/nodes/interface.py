import bpy
from mathutils import Vector
from .arrange import arrange_tree

NODE_SPACING = 250


def input_named_attribute(
    socket: bpy.types.NodeSocket, name: str, data_type: str | None
) -> bpy.types.GeometryNode:
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


def remove_linked(socket: bpy.types.NodeSocket) -> None:
    if socket.is_linked:
        socket.node.id_data.links.remove(socket.links[0])
        arrange_tree(socket.node.id_data)
