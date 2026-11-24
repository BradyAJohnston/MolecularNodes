from typing import Any, Iterable, List, Literal
import bpy
from bpy.types import GeometryNodeTree, Node
from nodebpy import TreeBuilder, arrange
from nodebpy import geometry as g
from .. import color
from . import geometry as mng

STYLE_NODE_MAPPING = {
    "spheres": mng.StyleSpheres,
    "cartoon": mng.StyleCartoon,
    "ribbon": mng.StyleRibbon,
    "surface": mng.StyleSurface,
    "sticks": mng.StyleSticks,
    "ball_and_stick": mng.StyleBallAndStick,
}

STYLE_LITERALS = Literal[
    "spheres", "cartoon", "ribbon", "surface", "sticks", "ball_and_stick"
]

# style keyword -> node group that implements it, used when swapping an existing
# style node for a different style. Keys match the enum identifiers in
# ui.style.STYLE_ITEMS and ui.ops.DENSITY_STYLE_ITEMS
styles_mapping = {
    **STYLE_NODE_MAPPING,
    "density_surface": mng.DensityStyleSurface,
    "density_iso_surface": mng.DensityStyleISOSurface,
    "density_wire": mng.DensityStyleWire,
}


def get_output(tree: bpy.types.NodeTree) -> Node:
    "The Group Output node of the tree"
    return tree.nodes[bpy.app.translations.pgettext_data("Group Output")]


def swap(node: Node, tree: GeometryNodeTree) -> None:
    """
    Swap out the node's node_tree, keeping its links by socket name and its material.

    Other inputs are reset to the defaults of the new tree.
    """
    if not isinstance(node, bpy.types.GeometryNodeGroup):
        raise TypeError(
            f"Can only swap the tree of a GeometryNodeGroup, not {type(node).__name__}"
        )
    parent = node.id_data
    assert isinstance(parent, bpy.types.NodeTree)

    # links are removed before the swap, so Blender doesn't reconnect them itself
    input_links: list[tuple[Any, str]] = []
    for socket in node.inputs:
        for link in socket.links or ():
            input_links.append((link.from_socket, socket.name))
            parent.links.remove(link)
    output_links: list[tuple[str, Any]] = []
    for socket in node.outputs:
        for link in socket.links or ():
            output_links.append((socket.name, link.to_socket))
            parent.links.remove(link)
    material = getattr(node.inputs.get("Material"), "default_value", None)

    # only change the label if it hasn't been customised away from the tree name
    if node.node_tree is not None and node.label == node.node_tree.name:
        node.label = tree.name
    node.node_tree = tree
    node.name = tree.name

    for from_socket, name in input_links:
        if name in node.inputs:
            parent.links.new(from_socket, node.inputs[name])
    for name, to_socket in output_links:
        if name in node.outputs:
            parent.links.new(node.outputs[name], to_socket)

    assert tree.interface is not None
    for item in tree.interface.items_tree:
        if (
            isinstance(item, bpy.types.NodeTreeInterfaceSocket)
            and item.in_out == "INPUT"
            and hasattr(item, "default_value")
        ):
            socket: Any = node.inputs[item.identifier]
            socket.default_value = item.default_value
    if material is not None and "Material" in node.inputs:
        material_socket: Any = node.inputs["Material"]
        material_socket.default_value = material


def get_star_node(obj: bpy.types.Object) -> bpy.types.Node:
    "The Starfile Instances node in the object's Molecular Nodes tree"
    modifier = obj.modifiers["Molecular Nodes"]
    assert isinstance(modifier, bpy.types.NodesModifier)
    tree = modifier.node_group
    assert tree is not None
    for node in tree.nodes:
        node_tree = getattr(node, "node_tree", None)
        if node_tree is not None and "Starfile Instances" in node_tree.name:
            return node
    raise ValueError(f"No Starfile Instances node found in tree for {obj.name}")


def _linked_node(socket: bpy.types.NodeSocket) -> Node:
    "Get the node of the first link into this socket"
    links = socket.links
    assert links is not None
    node = links[0].from_node
    assert node is not None
    return node


def _previous_node(node: Node) -> Node:
    "Get the node which is the first connection to the first input of this node"
    return _linked_node(node.inputs[0])


def _final_join(tree: GeometryNodeTree) -> Node | None:
    """
    Get the last JoinGeometry node in the tree, if there is one.
    """
    try:
        current = _previous_node(get_output(tree))
        while True:
            if current.bl_idname == "GeometryNodeGroupInput":
                return None
            if current.bl_idname == "GeometryNodeJoinGeometry":
                return current
            current = _previous_node(current)
    except IndexError:
        return None


def get_final_style_nodes(
    tree: bpy.types.GeometryNodeTree,
) -> List[bpy.types.GeometryNodeGroup]:
    """
    Get the final style nodes in the tree.
    """
    join = _final_join(tree)
    node = join if join is not None else get_output(tree)
    links = node.inputs[0].links
    assert links is not None

    # substring rather than prefix match, so density styles ("Density Style
    # Surface") are found too, matching the styles list in the UI
    return [
        style
        for link in reversed(links)
        if isinstance(style := link.from_node, bpy.types.GeometryNodeGroup)
        and "Style" in style.name
    ]


def remove_style_node(node: Node) -> None:
    """
    Remove a style node from its tree, along with the nodes linked into its inputs.
    """
    tree = node.id_data
    assert isinstance(tree, bpy.types.NodeTree)
    to_remove = [node] + [
        _linked_node(input) for input in node.inputs if input.is_linked
    ]
    for node_to_remove in to_remove:
        tree.nodes.remove(node_to_remove)
    arrange(tree)


def custom_boolean_iswitch(
    name: str,
    items: Iterable[str],
    attribute_name: str = "chain_id",
    offset: int = 0,
    prefix: str = "",
) -> TreeBuilder[GeometryNodeTree]:
    with g.tree(name) as tree:
        attr = g.NamedAttribute.integer(attribute_name)

        switch = g.IndexSwitch.boolean(
            index=attr if offset == 0 else attr + offset,
            items=[tree.inputs.boolean(prefix + x) for x in items],
        )

        switch >> tree.outputs.boolean("Selection")
        ~switch >> tree.outputs.boolean("Inverted")

    tree.tree.color_tag = "INPUT"

    return tree


def custom_color_iswitch(
    name: str,
    items: dict[str, tuple[float, float, float, float]] | Iterable[int | float | str],
    attribute_name: str = "chain_id",
    offset: int = 0,
) -> TreeBuilder[GeometryNodeTree]:
    with g.tree(name) as tree:
        attr = g.NamedAttribute.integer(attribute_name)

        # random_rgb() returns an RGBA array rather than a tuple
        colors: dict[Any, Any] = (
            items
            if isinstance(items, dict)
            else {str(key): color.random_rgb() for key in items}
        )

        switch = g.IndexSwitch.color(
            index=attr if offset == 0 else attr + offset,
            items=[tree.inputs.color(key, value) for key, value in colors.items()],
        )

        switch >> tree.outputs.color()

    tree.tree.color_tag = "INPUT"
    return tree
