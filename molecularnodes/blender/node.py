from abc import ABCMeta
from typing import Optional, Any, Union, List

import bpy
from mathutils import Vector


class Tree:
    def __init__(self, tree: bpy.types.GeometryNodeTree):
        self.tree = tree

    @property
    def nodes(self):
        return self.tree.nodes

    def remove(self, node: bpy.types.GeometryNode) -> None:
        self.tree.nodes.remove(node=node)

    def add(
        self, name: str, location: Optional[List[float]] = None
    ) -> bpy.types.GeometryNode:
        if not name.startswith("GeometryNode"):
            name = "GeometryNode" + name

        node = self.tree.nodes.new(name)

        if location is not None:
            node.location = location

        return node

    def link(
        self,
        link_from: Union[bpy.types.NodeSocket, bpy.types.GeometryNode],
        link_to: Union[bpy.types.NodeSocket, bpy.types.GeometryNode],
    ):
        """
        Link two nodes together in the tree

        If the given inputs are not sockets, get the first sockets and attempt to link them,
        otherwise create a link using the given socket
        """
        if not isinstance(link_from, bpy.types.NodeSocket):
            link_from = link_from.outputs[0]

        if not isinstance(link_to, bpy.types.NodeSocket):
            link_to = link_to.inputs[0]

        self.tree.links.new(link_from, link_to)


class Node:
    def __init__(self, node: bpy.types.GeometryNode):
        self.node = node
        self.offset_width = 200

    @property
    def tree(self) -> Tree:
        return Tree(self.node.id_data)

    def link_to(
        self,
        node: bpy.types.GeoemtryNode,
        index_self: Optional[Union[int, str]] = None,
        index_node: Optional[Union[int, str]] = None,
    ) -> None:
        if index_self is None:
            index_self = 0
        if index_node is None:
            index_node = 0
        self.tree.link(self.node[index_self], node[index_node])

    def add_and_link(self, name: str):
        loc = Vector(self.node.location) + Vector([0, self.offset_width])
        return self.tree.add(name=name, location=loc)


class GeometryNode(Node):
    def __init__(self, node):
        super().__init__(node)

    def instances_to_points(
        self,
        selection: Optional[Union[bpy.types.NodeSocket, bpy.types.GeometryNode]],
        position: Optional[Union[bpy.types.NodeSocket, bpy.types.GeometryNode]],
        radius: Union[bpy.types.NodeSocket, bpy.types.GeometryNode, float] = 0.05,
    ) -> Node:
        node = self.add_and_link("GeometryNodeInstancesToPoints")

        if isinstance(radius, float):
            node.inputs["Radius"].default_value = radius


def testing_this_funciton(self: int) -> bool:
    if isinstance(self, int):
        return True
