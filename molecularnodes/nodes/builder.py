import bpy
from bpy.types import GeometryNode, GeometryNodeTree, NodeSocket


def source_socket(node: GeometryNode | NodeSocket):
    if isinstance(node, GeometryNode):
        return node.outputs[0]
    elif isinstance(node, NodeSocket):
        return node


def target_socket(node: GeometryNode | NodeSocket):
    if isinstance(node, GeometryNode):
        return node.inputs[0]
    elif isinstance(node, NodeSocket):
        return node


class TreeBuilder:
    def __init__(self, tree: bpy.types.GeometryNodeTree):
        self.tree = tree

    def link(self, socket1, socket2):
        self.tree.links.new(socket1, socket2)


class NodeBuilder:
    def __init__(self, node: bpy.types.GeometryNode):
        self.node = node

    @property
    def tree(self):
        return TreeBuilder(self.node.id_data)  # type: ignore

    def link_to(self, target: GeometryNode | NodeSocket):
        if isinstance(target, GeometryNode):
            self.tree.link(self.node.outputs[0], target.inputs[0])
        elif isinstance(target, NodeSocket):
            self.tree.link(self.node.outputs[0], target)

    def link_from(self, source: GeometryNode | NodeSocket, input: NodeSocket | str):
        source = source_socket(source)
        if isinstance(input, str):
            input = self.node.inputs[input]
        self.tree.link(source, input)
