import bpy
import databpy


def get_links(socket: bpy.types.NodeSocket) -> bpy.types.NodeLinks:
    """
    Get the links between two sockets
    """
    links = socket.links
    if links is None:
        raise ValueError("No links found for socket")
    return links


def int_menu_switch_from_dict(
    dictionary: dict, name: str = "NewMenuSwitch"
) -> bpy.types.NodeTree:
    """
    Create a menu from a dictionary of options
    """
    tree = databpy.nodes.new_tree(name=name, geometry=False)

    node_m = tree.nodes.new("GeometryNodeMenuSwitch")
    node_m.data_type = "INT"

    starting_enum_length = len(node_m.enum_items)

    for i, (key, value) in enumerate(dictionary.items()):
        if i < starting_enum_length:
            node_m.enum_items[i].name = key
        else:
            node_m.enum_items.new(name=key)

        node_m.inputs[i + 1].default_value = value

    tree.interface.new_socket("Selector", socket_type="NodeSocketMenu")
    tree.interface.new_socket("Value", socket_type="NodeSocketInt", in_out="OUTPUT")
    tree.links.new(databpy.nodes.get_input(tree).outputs["Selector"], node_m.inputs[0])
    tree.links.new(
        node_m.outputs["Output"], databpy.nodes.get_output(tree).inputs["Value"]
    )
    return tree
