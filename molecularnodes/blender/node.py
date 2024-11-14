from abc import ABCMeta
from typing import Optional, Any, Union, List
from bpy.types import NodeSocketBool, GeometryNode, NodeSocket

from .. import bpyd

import bpy
from mathutils import Vector


class Tree:
    def __init__(self, tree: bpy.types.GeometryNodeTree | str):
        if isinstance(tree, str):
            self.tree = bpy.data.node_groups.new(name=tree, type="GeometryNodeTree")
        else:
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
        link_from: NodeSocket | GeometryNode,
        link_to: NodeSocket | GeometryNode,
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
    def __init__(self, node: GeometryNode):
        self.node = node
        self.offset_width = 200

    @property
    def tree(self) -> Tree:
        return Tree(self.node.id_data)

    def link_to(
        self,
        node: GeometryNode,
        index_self: int | str = None,
        index_node: int | str = None,
    ) -> None:
        if index_self is None:
            index_self = 0
        if index_node is None:
            index_node = 0
        self.tree.link(self.node[index_self], node[index_node])


class TemporaryNodeTree:
    def __init__(self) -> None:
        self.tree = Tree(
            py.data.node_groups.new("TemporaryNodeGroup", "GeometryNodeTree")
        )

    def __enter__(self):
        return self.tree

    def __exit__(self):
        bpy.data.node_groups.remove(self.tree.tree)
        del self.tree


with TemporaryNodeTree() as tree:
    names = [name for name in dir(bpy.types) if name.startswith("GeometryNode")]
    for name in names:
        try:
            node = tree.add(name)
            print(f"Node: {name}")
            for input in node.inputs:
                print(f"  Input: {input.name}")
            for output in node.outputs:
                print(f"  Output: {output.name}")
        except RuntimeError:
            pass


# tree = node_group()
# (
#     tree.mesh_grid(x=1, y=20)
#     .instance_on_points(instance=icosphere().geometry)
#     .scale_instances()
#     .set_position(position=position() - 10)
#     .output()
# )
