from typing import List, Sequence
import bpy
from bpy.types import Node
from mathutils import Vector
from . import nodes
from .arrange import arrange_tree
from .interface import (
    TreeInterface,
    input_named_attribute,
    socket,
)
from .material import (
    assign_material,
    getset_material,
)
from .nodes import (
    NODE_SPACING,
    final_join,
    insert_before,
    loc_between,
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

    if isinstance(color, str) and color.lower() in ["default", "common"]:
        node_cc = insert_before(node_sc.inputs["Color"], "Color Common")
        node_car: bpy.types.GeometryNodeGroup = insert_before(  # type: ignore
            node_cc.inputs["Carbon"], "Color Attribute Random"
        )

        return node_car

    if isinstance(color, str):
        input_named_attribute(node_sc.inputs["Color"], color, "FLOAT_COLOR")
    else:
        node_sc.inputs["Color"].default_value = color  # type: ignore
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
        node_af.inputs["Frames"].default_value = frames  # type: ignore
    elif isinstance(frames, str):
        node_af.inputs["Frames"].default_value = bpy.data.collections[frames]  # type: ignore
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


def add_style_branch(
    tree: bpy.types.GeometryNodeTree,
    style: str | bpy.types.GeometryNodeTree,
    color: str | None = None,
    selection: str | None = None,
    material: bpy.types.Material | str | None = None,
    frames: bpy.types.Collection | str | None = None,
) -> None:
    """
    Add a style branch to the tree.
    """
    _frame = True
    link = tree.links.new
    input = nodes.get_input(tree)
    node_join = final_join(tree)

    current_min_y = min(node.location[1] for node in tree.nodes)
    ypos = current_min_y - 200
    xpos = loc_between(input, node_join, 0.75)[0]

    if isinstance(style, str):
        style_name = nodes.styles_mapping[style]
        nodes.append(style_name)
    elif isinstance(style, bpy.types.GeometryNodeTree):
        style_name = style.name
    else:
        raise ValueError(
            f"Style must be a string or a GeometryNodeTree, not {type(style)=}"
        )

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

    for nodelink in node_join.inputs[0].links:  # type: ignore
        if nodelink.from_socket.node == nodes.get_input(tree):  # type: ignore
            tree.links.remove(nodelink)

    # Apply style modifications
    if material:
        assign_material(node_style, material)
    if selection:
        input_named_attribute(node_style.inputs["Selection"], selection, "BOOLEAN")
    if color:
        insert_set_color(node_style, color)
    if frames:
        insert_animate_frames(node_style, frames)

    arrange_tree(tree)


def get_final_style_nodes(
    tree: bpy.types.GeometryNodeTree,
) -> List[bpy.types.Node | None]:
    """
    Get the final style nodes in the tree.
    """
    links: bpy.types.NodeLinks = final_join(tree).inputs[0].links  # type: ignore
    return [
        link.from_socket.node  # type: ignore
        for link in reversed(links)
        if link.from_socket.node.name.startswith("Style")  # type: ignore
    ]


class GeometryNodeInterFace(TreeInterface):
    """
    Interface for the geometry nodes in the tree.
    """

    def __init__(self, node: bpy.types.Node) -> None:
        super().__init__()
        self.tree: bpy.types.NodeTree = node.id_data
        self._nodes = []

    def remove(self) -> None:
        """
        Cleanup when this instance is explicitly deleted.
        """
        for node in self._nodes:
            self.tree.nodes.remove(node)
        arrange_tree(self.tree)

    def _expose_options(self, node: bpy.types.Node) -> None:
        self._nodes.append(node)
        for input in node.inputs:
            if not hasattr(input, "default_value"):
                continue

            # Create property name from node and input names
            prop_name = (
                "_".join([node.name.split(".")[0], input.name])
                .lower()
                .replace(" ", "_")
                .removeprefix("style_")
                .removeprefix("ribbon_")
                .removeprefix("surface_")
                .removeprefix("ball_and_stick_")
                .removeprefix("cartoon_")
                .removeprefix("spheres_")
            )
            if isinstance(input, bpy.types.NodeSocketMaterial):
                prop = getset_material(input)
            else:
                prop = socket(input)
            self._register_property(prop_name)
            setattr(self.__class__, prop_name, prop)


def create_style_interface(node: Node, linked: bool = True) -> GeometryNodeInterFace:
    """
    Dynamically create a StyleInterface class with exposed options.
    """
    class_name = f"DynamicStyleInterface_{node.name.replace('Style ', '')}"

    # Create the dynamic class with a unique name
    DynamicInterface = type(class_name, (GeometryNodeInterFace,), {})

    # Pre-expose the options for the main node
    interface = DynamicInterface(node)
    interface._expose_options(node)

    if not linked:
        return interface

    # Pre-expose options for linked nodes
    for input in node.inputs:
        if input.is_linked:
            linked_node = input.links[0].from_socket.node  # type: ignore
            interface._expose_options(linked_node)

    return interface


def style_interfaces_from_tree(
    tree: bpy.types.GeometryNodeTree,
) -> list[GeometryNodeInterFace]:
    return [create_style_interface(node) for node in get_final_style_nodes(tree)]
