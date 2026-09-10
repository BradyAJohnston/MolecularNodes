import bpy
from .arrange import arrange_tree

NODE_SPACING = 250


def remove_linked(socket: bpy.types.NodeSocket) -> None:
    if socket.is_linked:
        socket.node.id_data.links.remove(socket.links[0])
        arrange_tree(socket.node.id_data)
