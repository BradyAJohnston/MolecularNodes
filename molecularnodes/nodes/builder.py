from typing import Iterable, TypeVar
import bpy
from mathutils import Vector, Quaternion
import numpy as np
from bpy.types import (
    GeometryNode,
    GeometryNodeTree,
    NodeSocket,
    NodeSocketMenu,
    NodeSocketRotation,
    Nodes,
    NodeGroupOutput,
    NodeGroupInput,
    NodeSocketBool,
    NodeSocketVector,
)

GEO_NODE_NAMES = (
    f"GeometryNode{name}"
    for name in (
        "SetPosition",
        "TransformGeometry",
        "GroupInput",
        "GroupOutput",
        "MeshToPoints",
        "PointsToVertices",
    )
)

# POSSIBLE_NODE_NAMES = "GeometryNode"
LINKABLE = TypeVar("LINKABLE", bound=GeometryNode | NodeSocket | "NodeBuilder")
TYPE_INPUT_VECTOR = TypeVar("TYPE_INPUT_VECTOR", bound=NodeSocketVector | Vector | 'NodeBuilder'| list[float] | tuple[float, float, float] | None)
TYPE_INPUT_ROTATION = TypeVar("TYPE_INPUT_ROTATION", bound=NodeSocketRotation | Quaternion | 'NodeBuilder'| list[float] | tuple[float, float, float, float] | None)
TYPE_INPUT_BOOLEAN = TypeVar("TYPE_INPUT_BOOLEAN", bound=NodeSocketBool | bool | 'NodeBuilder'| None)

def source_socket(node: GeometryNode | NodeSocket | "NodeBuilder"):
    if isinstance(node, GeometryNode):
        return node.outputs[0]
    elif isinstance(node, NodeSocket):
        return node
    elif isinstance(node, NodeBuilder):
        return node.node.outputs[0]
    else:
        raise TypeError(f"Unsupported type: {type(node)}")


def target_socket(node: GeometryNode | NodeSocket | "NodeBuilder"):
    if isinstance(node, GeometryNode):
        return node.inputs[0]
    elif isinstance(node, NodeSocket):
        return node
    elif isinstance(node, NodeBuilder):
        return node.node.inputs[0]
    else:
        raise TypeError(f"Unsupported type: {type(node)}")

class NodeAdder:
    def __init__(self, tree: 'TreeBuilder'):
        self.tree = tree

    def Position(self) -> 'NodeBuilder':
        return Position(self.tree)

    def Vector(self, value: Vector | list[float] | tuple[float, float, float] | None = None) -> 'NodeBuilder':
        return Vector(self.tree, value)

class TreeBuilder:
    just_added: GeometryNode | None = None

    def __init__(self, tree: GeometryNodeTree | str | None = None):
        if isinstance(tree, str):
            self.tree = bpy.data.node_groups.new(tree, 'GeometryNodeTree')
        elif tree is None:
            self.tree = bpy.data.node_groups.new('GeometryNodeTree', 'GeometryNodeTree')
        else:
            assert isinstance(tree, GeometryNodeTree)
            self.tree = tree

        self.node_adder = NodeAdder(self)

    @property
    def nodes(self) -> Nodes:
        return self.tree.nodes

    def input(self) -> GeometryNode:
        try:
            return self.tree.nodes['Group Input'] # type: ignore
        except KeyError:
            return self.tree.nodes.new('NodeGroupInput') # type: ignore

    def input_node(self):
        return NodeBuilder(self.input())

    def output(self) -> GeometryNode:
        try:
            return self.tree.nodes['Group Output'] # type: ignore
        except KeyError:
            return self.tree.nodes.new('NodeGroupOutput') # type: ignore

    def link(self, socket1: NodeSocket, socket2: NodeSocket):
        self.tree.links.new(socket1, socket2)

    def add(self, name: str) -> GeometryNode:
        assert name in GEO_NODE_NAMES
        self.just_added = self.tree.nodes.new(name) # type: ignore
        assert self.just_added is not None
        return self.just_added




class NodeBuilder:
    node: GeometryNode
    _tree: TreeBuilder
    name: str

    def __init__(self, tree: TreeBuilder | GeometryNode):
        if isinstance(tree, TreeBuilder):
            self.tree = tree
        elif isinstance(tree, GeometryNode):
            self.tree = TreeBuilder(tree.id_data)  # type: ignore
            self.node = tree
        else:
            raise TypeError(f"Expected TreeBuilder or GeometryNode, got {type(tree)}")

        if self.name is None:
            self.name = self.node.bl_idname
        else:
            self.tree.add(self.name)

    @property
    def tree(self):
        if self.node is None:
            raise ValueError

        return TreeBuilder(self.node.id_data)  # type: ignore

    @tree.setter
    def tree(self, value: TreeBuilder):
        self._tree = value

    @property
    def default_input(self) -> NodeSocket:
        return self.node.inputs[0]

    @property
    def default_output(self) -> NodeSocket:
        return self.node.outputs[0]

    def link(self, source: LINKABLE, target: LINKABLE):
        self.tree.link(source_socket(source), target_socket(target))

    def link_to(self, target: GeometryNode | NodeSocket | "NodeBuilder"):
        self.tree.link(self.default_output, target_socket(target))

    def link_from(self, source: LINKABLE, input: LINKABLE | str):
        if isinstance(input, str):
            self.link(source, self.node.inputs[input])
        else:
            self.link(source, input)

    def _establish_links(self, **kwargs):
        for name, value in kwargs.items():
            if isinstance(value, (int, float, list, tuple, np.ndarray)):
                self.inputs[name].default_value = value
            else:
                self.link_from(value, name)


    def set_position(
        self,
        selection: TYPE_INPUT_BOOLEAN = None,
        position: TYPE_INPUT_VECTOR = None,
        offset: TYPE_INPUT_VECTOR = None,
    ) -> "NodeBuilder":
        node = SetPosition(self.tree)
        self.link_from(self.node, node.node.inputs[0])
        node._establish_links(selection=selection, position=position, offset=offset)
        return node

    def transform_geometry(
        self,
        method: NodeSocketMenu | None = None,
        location: TYPE_INPUT_VECTOR = None,
        rotation: TYPE_INPUT_ROTATION = None,
        scale: TYPE_INPUT_VECTOR = None,
    ) -> "NodeBuilder":
        node = TransformGeometry(self.tree)
        self.link_from(self.node, node.node.inputs[0])
        node._establish_links(method=method, location=location, rotation=rotation, scale=scale)
        return node

    def output(self) -> None:
        self.link_from(self.node, self.tree.output().inputs[0])

class TransformGeometry(NodeBuilder):
    name = "TransformGeometry"

    def __init__(
        self,
        tree: TreeBuilder,
        method: NodeSocketMenu | None = None,
        location: TYPE_INPUT_VECTOR = None,
        rotation: TYPE_INPUT_ROTATION = None,
        scale: TYPE_INPUT_VECTOR = None,
    ):
        super().__init__(tree)
        self._establish_links(method=method, location=location, rotation=rotation, scale=scale)

class SetPosition(NodeBuilder):
    name= "SetPosition"

class Position(NodeBuilder):
    name = "Position"

    @property
    def position(self):
        return self. node.outputs["Position"]

class Vector(NodeBuilder):
    name = "FunctionNodeInputVector"

    def __init__(
        self,
        tree: TreeBuilder,
        value: Vector | list[float] | tuple[float, float, float] = (0.0, 0.0, 0.0)
    ):
        super().__init__(tree)
        self.value = value

    @property
    def value(self) -> Vector:
        return self.inputs[0].default_value

    @value.setter
    def value(self, value: Vector | list[float] | tuple[float, float, float]):
        self.inputs[0].default_value = value

class DefaultTree:
    def __init__(self):
        self.tree = tree = TreeBuilder()
        input = self.tree.tree.interface.new_socket('Geometry', in_out='INPUT', socket_type='NodeSocketGeometry')
        output = self.tree.tree.interface.new_socket('Geometry', in_out='OUTPUT', socket_type='NodeSocketGeometry')
        self.tree.link(
            self.tree.input().outputs[0],
            self.tree.output().inputs[0]
        )


    def __enter__(self):
        return self.tree

    def __exit__(self, exc_type, exc_value, traceback):
        pass



with DefaultTree() as tree:
    a = tree.node_adder
    new_tree = (
        tree
        .input_node()
        .set_position(position=a.Position(), offset=a.Vector())
        .transform_geometry(location=(0, 0, 1))
        .output()
    )
