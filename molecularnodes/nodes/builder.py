from __future__ import annotations
from typing import TYPE_CHECKING, ClassVar
import bpy
import numpy as np
from bpy.types import (
    GeometryNodeTree,
    Node,
    Nodes,
    NodeSocket,
)
from .arrange import arrange_tree

if TYPE_CHECKING:
    from molecularnodes.nodes.sockets import SocketBase

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


def normalize_name(name: str) -> str:
    """Convert 'Geometry' or 'My Socket' to 'geometry' or 'my_socket'."""
    return name.lower().replace(" ", "_")


def denormalize_name(attr_name: str) -> str:
    """Convert 'geometry' or 'my_socket' to 'Geometry' or 'My Socket'."""
    return attr_name.replace("_", " ").title()


def source_socket(node: LINKABLE) -> NodeSocket:
    if isinstance(node, NodeSocket):
        return node
    elif isinstance(node, Node):
        return node.outputs[0]
    elif hasattr(node, "default_output"):
        # NodeBuilder or SocketNodeBuilder
        return node.default_output
    else:
        raise TypeError(f"Unsupported type: {type(node)}")


def target_socket(node: LINKABLE) -> NodeSocket:
    if isinstance(node, NodeSocket):
        return node
    elif isinstance(node, Node):
        return node.inputs[0]
    elif hasattr(node, "default_input"):
        # NodeBuilder or SocketNodeBuilder
        return node.default_input
    else:
        raise TypeError(f"Unsupported type: {type(node)}")


class TreeBuilder:
    """Builder for creating Blender geometry node trees with a clean Python API."""

    _active_tree: ClassVar["TreeBuilder | None"] = None
    just_added: "Node | None" = None

    def __init__(self, tree: "GeometryNodeTree | str | None" = None):
        if isinstance(tree, str):
            self.tree = bpy.data.node_groups.new(tree, "GeometryNodeTree")
        elif tree is None:
            self.tree = bpy.data.node_groups.new("GeometryNodeTree", "GeometryNodeTree")
        else:
            assert isinstance(tree, GeometryNodeTree)
            self.tree = tree

        # Create socket accessors for named access
        self.inputs = SocketAccessor(self, "INPUT")
        self.outputs = SocketAccessor(self, "OUTPUT")

    def __enter__(self):
        TreeBuilder._active_tree = self
        return self

    def __exit__(self, *args):
        self.arrange()
        TreeBuilder._active_tree = None

    @property
    def nodes(self) -> Nodes:
        return self.tree.nodes

    def arrange(self):
        arrange_tree(self.tree)

    def _input_node(self) -> Node:
        """Get or create the Group Input node."""
        try:
            return self.tree.nodes["Group Input"]  # type: ignore
        except KeyError:
            return self.tree.nodes.new("NodeGroupInput")  # type: ignore

    def _output_node(self) -> Node:
        """Get or create the Group Output node."""
        try:
            return self.tree.nodes["Group Output"]  # type: ignore
        except KeyError:
            return self.tree.nodes.new("NodeGroupOutput")  # type: ignore

    def link(self, socket1: NodeSocket, socket2: NodeSocket):
        self.tree.links.new(socket1, socket2)

    def add(self, name: str) -> Node:
        self.just_added = self.tree.nodes.new(name)  # type: ignore
        assert self.just_added is not None
        return self.just_added

    def interface(
        self,
        inputs: list[SocketBase] | None = None,
        outputs: list[SocketBase] | None = None,
    ) -> None:
        """Define the node group interface with typed socket definitions.

        Args:
            inputs: List of input socket definitions
            outputs: List of output socket definitions

        Example:
            from molecularnodes.nodes.sockets import SocketGeometry, SocketBoolean, SocketVector

            tree.interface(
                inputs=[
                    SocketGeometry(name="Geometry"),
                    SocketBoolean(name="Selection", default=True),
                    SocketVector(name="Offset", default=(0, 0, 0)),
                ],
                outputs=[
                    SocketGeometry(name="Geometry"),
                ]
            )
        """
        if inputs:
            for socket_def in inputs:
                self._create_input_socket(socket_def)

        if outputs:
            for socket_def in outputs:
                self._create_output_socket(socket_def)

    def _create_input_socket(self, socket_def: SocketBase) -> None:
        """Create an input socket from a socket definition."""
        socket = self.tree.interface.new_socket(
            name=socket_def.name,
            in_out="INPUT",
            socket_type=socket_def.bl_socket_type,
        )
        self._configure_socket(socket, socket_def)

    def _create_output_socket(self, socket_def: SocketBase) -> None:
        """Create an output socket from a socket definition."""
        socket = self.tree.interface.new_socket(
            name=socket_def.name,
            in_out="OUTPUT",
            socket_type=socket_def.bl_socket_type,
        )
        self._configure_socket(socket, socket_def)

    def _configure_socket(self, socket, socket_def: SocketBase) -> None:
        """Configure socket properties from definition."""
        # Set default value if it exists
        if hasattr(socket_def, "default") and hasattr(socket, "default_value"):
            socket.default_value = socket_def.default

        # Set min/max values if they exist
        if hasattr(socket_def, "min_value") and socket_def.min_value is not None:
            if hasattr(socket, "min_value"):
                socket.min_value = socket_def.min_value

        if hasattr(socket_def, "max_value") and socket_def.max_value is not None:
            if hasattr(socket, "max_value"):
                socket.max_value = socket_def.max_value

        # Set description
        if socket_def.description:
            socket.description = socket_def.description


class NodeBuilder:
    """Base class for all geometry node wrappers."""

    node: Node
    _tree: "TreeBuilder"
    name: str
    _link_target: str | None = None  # Track which input should receive links

    def __init__(self):
        # Get active tree from context manager
        tree = TreeBuilder._active_tree
        if tree is None:
            raise RuntimeError(
                f"Node '{self.__class__.__name__}' must be created within a TreeBuilder context manager.\n"
                f"Usage:\n"
                f"  with tree:\n"
                f"      node = {self.__class__.__name__}()\n"
            )

        self._tree = tree
        self._link_target = None
        if hasattr(self.__class__, "name") and self.__class__.name is not None:
            self.node = self._tree.add(self.__class__.name)
        else:
            raise ValueError(
                f"Class {self.__class__.__name__} must define a 'name' attribute"
            )

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
            if value is ...:
                # Ellipsis indicates this input should receive links from >> operator
                self._link_target = name
                continue
            if isinstance(value, (int, float, list, tuple, np.ndarray)):
                try:
                    self.node.inputs[name].default_value = value
                except KeyError:
                    input = name.replace("_", " ").title()
                    self.node.inputs[input].default_value = value
            else:
                self.link_from(value, name)

    def __rshift__(self, other: "NodeBuilder") -> "NodeBuilder":
        """Chain nodes using >> operator. Links output to input.

        Usage:
            node1 >> node2 >> node3
            tree.inputs.value >> Math.add_(..., 0.1) >> tree.outputs.result

        If the target node has an ellipsis placeholder (...), links to that specific input.
        Otherwise, tries to find Geometry sockets first, then falls back to default.

        Returns the right-hand node to enable continued chaining.
        """
        # Get source socket
        try:
            self_out = self.node.outputs.get("Geometry") or self.default_output
        except (KeyError, IndexError):
            self_out = self.default_output

        # Get target socket - use link target if specified by ellipsis
        if other._link_target is not None:
            try:
                other_in = other.node.inputs[other._link_target]
            except KeyError:
                # Try with title case if direct access fails
                target_name = other._link_target.replace("_", " ").title()
                other_in = other.node.inputs[target_name]
        else:
            # Default behavior - try Geometry first, then default input
            try:
                other_in = other.node.inputs.get("Geometry") or other.default_input
            except (KeyError, IndexError):
                other_in = other.default_input

        self.tree.link(self_out, other_in)
        return other


class SocketNodeBuilder(NodeBuilder):
    """Special NodeBuilder for accessing specific sockets on input/output nodes."""

    def __init__(self, node: Node, socket_name: str, direction: str):
        # Don't call super().__init__ - we already have a node
        self.node = node
        self._tree = TreeBuilder(node.id_data)  # type: ignore
        self._socket_name = socket_name
        self._direction = direction

    @property
    def default_output(self) -> NodeSocket:
        """Return the specific named output socket."""
        if self._direction == "INPUT":
            return self.node.outputs[self._socket_name]
        else:
            # Output nodes don't have outputs, this shouldn't be called
            return self.node.outputs[0]

    @property
    def default_input(self) -> NodeSocket:
        """Return the specific named input socket."""
        if self._direction == "OUTPUT":
            return self.node.inputs[self._socket_name]
        else:
            # Input nodes don't have inputs, this shouldn't be called
            return self.node.inputs[0]


class SocketAccessor:
    """Provides named access to tree input or output sockets.

    Usage:
        tree.inputs.geometry  # Access "Geometry" input socket
        tree.outputs.result   # Access "Result" output socket
    """

    def __init__(self, tree: "TreeBuilder", direction: str):
        self._tree = tree
        self._direction = direction  # 'INPUT' or 'OUTPUT'

    def __getattr__(self, name: str) -> "NodeBuilder":
        """Access a socket by normalized name.

        Example:
            tree.inputs.geometry -> accesses "Geometry" socket
            tree.inputs.my_socket -> accesses "My Socket" socket
        """
        # Convert attribute name to socket name
        socket_name = denormalize_name(name)

        # Get the appropriate node
        if self._direction == "INPUT":
            node = self._tree._input_node()
        else:
            node = self._tree._output_node()

        # Return a NodeBuilder wrapping this specific socket access
        # We create a wrapper that will connect the right socket
        return SocketNodeBuilder(node, socket_name, self._direction)
