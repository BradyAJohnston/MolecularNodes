from __future__ import annotations
from typing import Iterable
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
    Node,
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
LINKABLE = "Node | NodeSocket | NodeBuilder"
TYPE_INPUT_VECTOR = "NodeSocketVector | Vector | NodeBuilder | list[float] | tuple[float, float, float] | None"
TYPE_INPUT_ROTATION = "NodeSocketRotation | Quaternion | NodeBuilder | list[float] | tuple[float, float, float, float] | None"
TYPE_INPUT_BOOLEAN = "NodeSocketBool | bool | NodeBuilder | None"

def source_socket(node: LINKABLE) -> NodeSocket:
    if isinstance(node, Node):
        return node.outputs[0]
    elif isinstance(node, NodeSocket):
        return node
    elif isinstance(node, NodeBuilder):
        return node.node.outputs[0]
    else:
        raise TypeError(f"Unsupported type: {type(node)}")


def target_socket(node: LINKABLE) -> NodeSocket:
    if isinstance(node, Node):
        return node.inputs[0]
    elif isinstance(node, NodeSocket):
        return node
    elif isinstance(node, NodeBuilder):
        return node.node.inputs[0]
    else:
        raise TypeError(f"Unsupported type: {type(node)}")

class NodeAdder:
    def __init__(self, tree: "TreeBuilder"):
        self.tree = tree

    def Position(self) -> "NodeBuilder":
        return Position(self.tree)

    def Vector(self, value: "Vector | list[float] | tuple[float, float, float] | None" = None) -> "NodeBuilder":
        return Vector(self.tree, value)

class TreeBuilder:
    just_added: "Node | None" = None

    def __init__(self, tree: "GeometryNodeTree | str | None" = None):
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

    def input(self) -> Node:
        try:
            return self.tree.nodes['Group Input'] # type: ignore
        except KeyError:
            return self.tree.nodes.new('NodeGroupInput') # type: ignore

    def input_node(self):
        return NodeBuilder(self.input())

    def output(self) -> Node:
        try:
            return self.tree.nodes['Group Output'] # type: ignore
        except KeyError:
            return self.tree.nodes.new('NodeGroupOutput') # type: ignore

    def link(self, socket1: NodeSocket, socket2: NodeSocket):
        self.tree.links.new(socket1, socket2)

    def add(self, name: str) -> Node:
        self.just_added = self.tree.nodes.new(name) # type: ignore
        assert self.just_added is not None
        return self.just_added




class NodeBuilder:
    node: Node
    _tree: "TreeBuilder"
    name: str

    def __init__(self, tree: "TreeBuilder | Node"):
        if isinstance(tree, TreeBuilder):
            self._tree = tree
            if hasattr(self.__class__, 'name') and self.__class__.name is not None:
                self.node = self._tree.add(self.__class__.name)
            else:
                raise ValueError(f"Class {self.__class__.__name__} must define a 'name' attribute")
        elif isinstance(tree, Node):  # Check if it's a Blender node type
            self._tree = TreeBuilder(tree.id_data)  # type: ignore
            self.node = tree
        else:
            raise TypeError(f"Expected TreeBuilder or Node, got {type(tree)}")

    @property
    def tree(self) -> "TreeBuilder":
        return self._tree

    @tree.setter
    def tree(self, value: "TreeBuilder"):
        self._tree = value

    @property
    def default_input(self) -> NodeSocket:
        return self.node.inputs[0]

    @property
    def default_output(self) -> NodeSocket:
        return self.node.outputs[0]

    def link(self, source: LINKABLE, target: LINKABLE):
        self.tree.link(source_socket(source), target_socket(target))

    def link_to(self, target: LINKABLE):
        self.tree.link(self.default_output, target_socket(target))

    def link_from(self, source: LINKABLE, input: "LINKABLE | str"):
        if isinstance(input, str):
            try:
                self.link(source, self.node.inputs[input])
            except KeyError:
                input = input.replace("_", " ").title()
                self.link(source, self.node.inputs[input])
        else:
            self.link(source, input)

    def _establish_links(self, **kwargs):
        for name, value in kwargs.items():
            if value is None:
                continue
            if isinstance(value, (int, float, list, tuple, np.ndarray)):
                try:
                    self.node.inputs[name].default_value = value
                except KeyError:
                    input = name.replace("_", " ").title()
                    self.node.inputs[input].default_value = value
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
        method: "NodeSocketMenu | None" = None,
        translation: TYPE_INPUT_VECTOR = None,
        rotation: TYPE_INPUT_ROTATION = None,
        scale: TYPE_INPUT_VECTOR = None,
    ) -> "NodeBuilder":
        node = TransformGeometry(self.tree)
        self.link_from(self.node, node.node.inputs[0])
        node._establish_links(method=method, translation=translation, rotation=rotation, scale=scale)
        return node

    def output(self) -> None:
        self.link_from(self.node, self.tree.output().inputs[0])

class TransformGeometry(NodeBuilder):
    name = "GeometryNodeTransform"

    def __init__(
        self,
        tree: "TreeBuilder",
        method: "NodeSocketMenu | None" = None,
        translation: TYPE_INPUT_VECTOR = None,
        rotation: TYPE_INPUT_ROTATION = None,
        scale: TYPE_INPUT_VECTOR = None,
    ):
        super().__init__(tree)
        self._establish_links(method=method, translation=translation, rotation=rotation, scale=scale)

class SetPosition(NodeBuilder):
    name = "GeometryNodeSetPosition"

class Position(NodeBuilder):
    name = "GeometryNodeInputPosition"

    @property
    def position(self):
        return self.node.outputs["Position"]

class Vector(NodeBuilder):
    name = "FunctionNodeInputVector"

    def __init__(
        self,
        tree: "TreeBuilder",
        value: "Vector | list[float] | tuple[float, float, float]" = (0.0, 0.0, 0.0)
    ):
        super().__init__(tree)
        if value is None:
            return
        self.value = value

    @property
    def value(self) -> Vector:
        node: bpy.types.FunctionNodeInputVector = self.node
        return node.vector

    @value.setter
    def value(self, value: "Vector | list[float] | tuple[float, float, float]"):
        node: bpy.types.FunctionNodeInputVector = self.node
        node.vector = value

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
        .transform_geometry(translation=(0, 0, 1))
        .output()
    )
