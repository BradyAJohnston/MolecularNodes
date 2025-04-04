def get_links(socket: bpy.types.NodeSocket) -> bpy.types.NodeLinks:
    """
    Get the links between two sockets
    """
    links = socket.links
    if links is None:
        raise ValueError("No links found for socket")
    return links
