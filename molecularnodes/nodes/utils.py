from typing import Iterable, List, Literal
import bpy
from bpy.types import GeometryNodeTree, Node
from databpy.nodes import append_from_blend, get_output, swap_tree
from nodebpy import TreeBuilder
from nodebpy import geometry as g
from nodebpy.builder.arrange import arrange_tree
from .. import color
from ..assets import MN_DATA_FILE
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

# style keyword -> name of the node tree that implements it, used when swapping
# an existing style node for a different style. Keys match the enum identifiers
# in ui.style.STYLE_ITEMS and ui.ops.DENSITY_STYLE_ITEMS
styles_mapping = {
    "preset_1": "Style Preset 1",
    "preset_2": "Style Preset 2",
    "preset_3": "Style Preset 3",
    "preset_4": "Style Preset 4",
    "spheres": "Style Spheres",
    "cartoon": "Style Cartoon",
    "sticks": "Style Sticks",
    "ribbon": "Style Ribbon",
    "surface": "Style Surface",
    "ball_and_stick": "Style Ball and Stick",
    "density_surface": "Density Style Surface",
    "density_iso_surface": "Density Style ISO Surface",
    "density_wire": "Density Style Wire",
}


def append(name: str, link: bool = False) -> bpy.types.GeometryNodeTree:
    "Append a GN node group from the MN data file"
    GN_TREES_PATH = MN_DATA_FILE / "NodeTree"
    return append_from_blend(name, filepath=str(GN_TREES_PATH), link=link)


def swap(node: bpy.types.Node, tree: str | bpy.types.NodeTree) -> None:
    "Swap out the node's node_tree, maintaining the old socket connections"
    if isinstance(tree, str):
        try:
            tree = bpy.data.node_groups[tree]
        except KeyError:
            tree = append(tree)
    # only change the label if it hasn't been customised away from the tree name
    if node.label == node.node_tree.name:
        node.label = tree.name
    swap_tree(node=node, tree=tree)


def get_star_node(obj: bpy.types.Object) -> bpy.types.Node:
    "The Starfile Instances node in the object's Molecular Nodes tree"
    tree = obj.modifiers["Molecular Nodes"].node_group
    for node in tree.nodes:
        node_tree = getattr(node, "node_tree", None)
        if node_tree is not None and "Starfile Instances" in node_tree.name:
            return node
    raise ValueError(f"No Starfile Instances node found in tree for {obj.name}")


def _previous_node(node: Node) -> Node:
    "Get the node which is the first connection to the first input of this node"
    return node.inputs[0].links[0].from_socket.node


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
    links: bpy.types.NodeLinks = node.inputs[0].links

    # substring rather than prefix match, so density styles ("Density Style
    # Surface") are found too, matching the styles list in the UI
    return [
        link.from_socket.node
        for link in reversed(links)
        if "Style" in link.from_socket.node.name
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

        if not isinstance(items, dict):
            items: dict[str, tuple[float, ...]] = {
                str(key): color.random_rgb() for key in items
            }

        switch = g.IndexSwitch.color(
            index=attr if offset == 0 else attr + offset,
            items=[tree.inputs.color(key, value) for key, value in items.items()],
        )

        switch >> tree.outputs.color()

    tree.tree.color_tag = "INPUT"
    return tree
