"""
Node tree builder providing fluent API for creating Blender node trees.
"""

from typing import Any, Literal, Union
from mathutils import Vector

import bpy

from .node import NodeWrapper, GeometryNodeGroupWrapper
from .socket import SocketWrapper, SocketCollection


# Socket type mapping
SOCKET_TYPES = {
    "BOOLEAN": "NodeSocketBool",
    "GEOMETRY": "NodeSocketGeometry",
    "INT": "NodeSocketInt",
    "MATERIAL": "NodeSocketMaterial",
    "VECTOR": "NodeSocketVector",
    "STRING": "NodeSocketString",
    "VALUE": "NodeSocketFloat",
    "FLOAT": "NodeSocketFloat",
    "COLLECTION": "NodeSocketCollection",
    "TEXTURE": "NodeSocketTexture",
    "COLOR": "NodeSocketColor",
    "RGBA": "NodeSocketColor",
    "IMAGE": "NodeSocketImage",
}


class NodeTreeBuilder:
    """
    Fluent builder for Blender node trees.

    Examples:
        # Create geometry node tree
        tree = NodeTreeBuilder.geometry("My Nodes")

        # Add nodes
        compare = tree.add_node("FunctionNodeCompare",
            name="Compare X",
            operation="LESS_EQUAL",
            location=(100, 200)
        )

        # Link nodes
        tree.link(output_socket, input_socket)

        # Get the final tree
        bpy_tree = tree.build()
    """

    def __init__(self, name: str, tree_type: str = "GeometryNodeTree"):
        """
        Initialize node tree builder.

        Args:
            name: Name for the node tree
            tree_type: Type of node tree ("GeometryNodeTree", "ShaderNodeTree", "CompositorNodeTree")
        """
        self._name = name
        self._tree_type = tree_type
        self._tree: bpy.types.NodeTree | None = None
        self._input_node: NodeWrapper | None = None
        self._output_node: NodeWrapper | None = None
        self._fallback = True

    @classmethod
    def create(cls, name: str) -> "NodeTreeBuilder":
        """
        Create a new node tree builder.

        Args:
            name: Name for the node tree

        Returns:
            New NodeTreeBuilder instance
        """
        return cls(name)

    @classmethod
    def geometry(
        cls,
        name: str = "Geometry Nodes",
        fallback: bool = True,
        input_name: str = "Geometry",
        output_name: str = "Geometry",
    ) -> "NodeTreeBuilder":
        """
        Create a geometry node tree with standard setup.

        Args:
            name: Name for the node tree
            fallback: If True, return existing tree with same name if it exists
            input_name: Name for the input geometry socket
            output_name: Name for the output geometry socket

        Returns:
            Configured NodeTreeBuilder instance
        """
        builder = cls(name, "GeometryNodeTree")
        builder._fallback = fallback
        builder._init_geometry_tree(input_name, output_name)
        return builder

    @classmethod
    def shader(cls, name: str = "Shader Nodes") -> "NodeTreeBuilder":
        """
        Create a shader node tree.

        Args:
            name: Name for the node tree

        Returns:
            Configured NodeTreeBuilder instance
        """
        builder = cls(name, "ShaderNodeTree")
        builder._init_tree()
        return builder

    @classmethod
    def compositor(cls, name: str = "Compositor Nodes") -> "NodeTreeBuilder":
        """
        Create a compositor node tree.

        Args:
            name: Name for the node tree

        Returns:
            Configured NodeTreeBuilder instance
        """
        builder = cls(name, "CompositorNodeTree")
        builder._init_tree()
        return builder

    def _init_tree(self) -> "NodeTreeBuilder":
        """Initialize the underlying Blender node tree."""
        # Check if tree exists
        if self._fallback:
            existing = bpy.data.node_groups.get(self._name)
            if existing:
                self._tree = existing
                return self

        # Create new tree
        self._tree = bpy.data.node_groups.new(name=self._name, type=self._tree_type)
        return self

    def _init_geometry_tree(
        self, input_name: str = "Geometry", output_name: str = "Geometry"
    ) -> "NodeTreeBuilder":
        """
        Initialize geometry node tree with input/output nodes and sockets.

        Args:
            input_name: Name for input geometry socket
            output_name: Name for output geometry socket

        Returns:
            Self for chaining
        """
        self._init_tree()

        # Create input/output nodes
        input_node = self._tree.nodes.new("NodeGroupInput")
        output_node = self._tree.nodes.new("NodeGroupOutput")

        # Position them
        input_node.location.x = -200 - input_node.width
        output_node.location.x = 200

        # Create geometry sockets
        self._tree.interface.new_socket(
            input_name, in_out="INPUT", socket_type="NodeSocketGeometry"
        )
        self._tree.interface.new_socket(
            output_name, in_out="OUTPUT", socket_type="NodeSocketGeometry"
        )

        # Link input to output by default
        self._tree.links.new(output_node.inputs[0], input_node.outputs[0])

        # Wrap nodes
        self._input_node = NodeWrapper(self, "NodeGroupInput", input_node)
        self._output_node = NodeWrapper(self, "NodeGroupOutput", output_node)

        return self

    @property
    def tree(self) -> bpy.types.NodeTree:
        """Get the underlying Blender node tree."""
        if self._tree is None:
            self._init_tree()
        return self._tree

    @property
    def name(self) -> str:
        """Get the tree name."""
        return self._tree.name if self._tree else self._name

    @name.setter
    def name(self, value: str):
        """Set the tree name."""
        if self._tree:
            self._tree.name = value
        self._name = value

    @property
    def input_node(self) -> NodeWrapper | None:
        """Get the Group Input node (if it exists)."""
        return self._input_node

    @property
    def output_node(self) -> NodeWrapper | None:
        """Get the Group Output node (if it exists)."""
        return self._output_node

    @property
    def inputs(self) -> SocketCollection | None:
        """Get the input sockets of the Group Input node."""
        if self._input_node:
            return self._input_node.outputs  # Group Input outputs are the tree inputs
        return None

    @property
    def outputs(self) -> SocketCollection | None:
        """Get the output sockets of the Group Output node."""
        if self._output_node:
            return self._output_node.inputs  # Group Output inputs are the tree outputs
        return None

    def add_input(
        self,
        name: str,
        socket_type: str,
        description: str = "",
        default_value: Any = None,
        **kwargs,
    ) -> "NodeTreeBuilder":
        """
        Add an input socket to the tree interface.

        Args:
            name: Socket name
            socket_type: Socket type (e.g., "GEOMETRY", "FLOAT", "INT")
            description: Socket description
            default_value: Default value for the socket
            **kwargs: Additional socket properties (min_value, max_value, etc.)

        Returns:
            Self for chaining
        """
        socket_type_name = SOCKET_TYPES.get(socket_type.upper(), socket_type)

        socket = self._tree.interface.new_socket(
            name, in_out="INPUT", socket_type=socket_type_name
        )

        if description:
            socket.description = description

        if default_value is not None:
            socket.default_value = default_value

        # Set additional properties
        for key, value in kwargs.items():
            setattr(socket, key, value)

        return self

    def add_output(
        self,
        name: str,
        socket_type: str,
        description: str = "",
        **kwargs,
    ) -> "NodeTreeBuilder":
        """
        Add an output socket to the tree interface.

        Args:
            name: Socket name
            socket_type: Socket type (e.g., "GEOMETRY", "FLOAT", "INT")
            description: Socket description
            **kwargs: Additional socket properties

        Returns:
            Self for chaining
        """
        socket_type_name = SOCKET_TYPES.get(socket_type.upper(), socket_type)

        socket = self._tree.interface.new_socket(
            name, in_out="OUTPUT", socket_type=socket_type_name
        )

        if description:
            socket.description = description

        # Set additional properties
        for key, value in kwargs.items():
            setattr(socket, key, value)

        return self

    def add_inputs(self, inputs: dict[str, tuple | Any]) -> "NodeTreeBuilder":
        """
        Add multiple input sockets at once.

        Args:
            inputs: Dict mapping socket name to (type, config) or (type, default_value)

        Returns:
            Self for chaining

        Example:
            tree.add_inputs({
                "Geometry": ("GEOMETRY", "Input geometry"),
                "Threshold": ("FLOAT", {"default": 0.5, "min": 0, "max": 1}),
                "Visible": ("BOOL", True),
            })
        """
        for name, config in inputs.items():
            if isinstance(config, tuple):
                socket_type = config[0]
                if len(config) > 1:
                    extra = config[1]
                    if isinstance(extra, dict):
                        # Dict with properties
                        self.add_input(name, socket_type, **extra)
                    else:
                        # Simple default value or description
                        if isinstance(extra, str):
                            self.add_input(name, socket_type, description=extra)
                        else:
                            self.add_input(name, socket_type, default_value=extra)
                else:
                    self.add_input(name, socket_type)
            else:
                # Just a type string
                self.add_input(name, config)

        return self

    def add_node(
        self,
        node_type: str,
        name: str | None = None,
        location: tuple[float, float] | None = None,
        **properties,
    ) -> NodeWrapper:
        """
        Add a node to the tree.

        Args:
            node_type: Blender node type (e.g., "FunctionNodeCompare")
            name: Optional node name
            location: Optional (x, y) location
            **properties: Node properties to set (operation, data_type, etc.)

        Returns:
            NodeWrapper for the created node

        Example:
            compare = tree.add_node("FunctionNodeCompare",
                name="Compare X",
                operation="LESS_EQUAL",
                location=(100, 200)
            )
        """
        node = NodeWrapper(self, node_type)

        if name:
            node.set_name(name)

        if location:
            node.set_location(*location)

        if properties:
            node.set_properties(**properties)

        return node

    def add_group(
        self,
        group_name: str,
        name: str | None = None,
        location: tuple[float, float] | None = None,
        material: str | bpy.types.Material | None = None,
        link: bool = False,
        subtree_name: str | None = None,
        **properties,
    ) -> GeometryNodeGroupWrapper:
        """
        Add a GeometryNodeGroup (custom node group) to the tree.

        Args:
            group_name: Name of the node group to load (from MN data file)
            name: Optional node name (defaults to group_name)
            location: Optional (x, y) location
            material: Optional material to assign
            link: Whether to link the node group
            subtree_name: Optional name for the node tree instance
            **properties: Additional node properties

        Returns:
            GeometryNodeGroupWrapper for the created node

        Example:
            style = tree.add_group("Style Density Surface",
                location=(400, 0),
                material="default"
            )
        """
        # Import here to avoid circular dependency
        from ...nodes.nodes import append

        node_group = GeometryNodeGroupWrapper(self)

        # Load the node tree
        node_group.node_tree = append(group_name, link=link)

        if name:
            node_group.set_name(name)
        else:
            node_group.set_name(group_name)

        if location:
            node_group.set_location(*location)

        if material:
            node_group.with_material(material)

        if subtree_name:
            node_group.set_subtree_name(subtree_name)

        if properties:
            node_group.set_properties(**properties)

        return node_group

    def link(
        self,
        from_socket: Union[SocketWrapper, bpy.types.NodeSocket],
        to_socket: Union[SocketWrapper, bpy.types.NodeSocket],
    ) -> "NodeTreeBuilder":
        """
        Create a link between two sockets.

        Args:
            from_socket: Output socket (source)
            to_socket: Input socket (destination)

        Returns:
            Self for chaining
        """
        # Unwrap sockets if needed
        from_sock = from_socket.socket if isinstance(from_socket, SocketWrapper) else from_socket
        to_sock = to_socket.socket if isinstance(to_socket, SocketWrapper) else to_socket

        self._tree.links.new(to_sock, from_sock)
        return self

    def connect_chain(self, nodes: list[NodeWrapper | bpy.types.Node]) -> "NodeTreeBuilder":
        """
        Connect a chain of nodes together (output[0] -> input[0]).

        Args:
            nodes: List of nodes to connect in sequence

        Returns:
            Self for chaining

        Example:
            tree.connect_chain([input_node, process_node, join_node, output_node])
        """
        for i in range(len(nodes) - 1):
            from_node = nodes[i]
            to_node = nodes[i + 1]

            # Unwrap if needed
            from_n = from_node.node if isinstance(from_node, NodeWrapper) else from_node
            to_n = to_node.node if isinstance(to_node, NodeWrapper) else to_node

            # Link first output to first input
            self._tree.links.new(to_n.inputs[0], from_n.outputs[0])

        return self

    def auto_layout(self, spacing: int = 50, add_group_input: bool = True) -> "NodeTreeBuilder":
        """
        Automatically arrange nodes in the tree.

        Args:
            spacing: Spacing between nodes
            add_group_input: Whether to include group input in layout

        Returns:
            Self for chaining
        """
        # Import here to avoid circular dependency
        from ...nodes.arrange import arrange_tree

        arrange_tree(self._tree, spacing=(spacing, 25), add_group_input=add_group_input)
        return self

    def get_node(self, name: str) -> NodeWrapper:
        """
        Get a node by name.

        Args:
            name: Node name

        Returns:
            Wrapped node

        Raises:
            KeyError: If node not found
        """
        node = self._tree.nodes.get(name)
        if node is None:
            raise KeyError(f"Node '{name}' not found in tree")

        # Wrap appropriately
        if isinstance(node, bpy.types.GeometryNodeGroup):
            return GeometryNodeGroupWrapper(self, node)
        else:
            return NodeWrapper(self, node.bl_idname, node)

    def remove_node(self, node: NodeWrapper | bpy.types.Node | str) -> "NodeTreeBuilder":
        """
        Remove a node from the tree.

        Args:
            node: Node to remove (wrapper, raw node, or name)

        Returns:
            Self for chaining
        """
        if isinstance(node, str):
            node_obj = self._tree.nodes[node]
        elif isinstance(node, NodeWrapper):
            node_obj = node.node
        else:
            node_obj = node

        self._tree.nodes.remove(node_obj)
        return self

    def clear_nodes(self) -> "NodeTreeBuilder":
        """
        Remove all nodes from the tree.

        Returns:
            Self for chaining
        """
        self._tree.nodes.clear()
        return self

    def build(self) -> bpy.types.NodeTree:
        """
        Finalize and return the underlying Blender node tree.

        Returns:
            The bpy.types.NodeTree object
        """
        return self._tree

    def __repr__(self) -> str:
        return f"<NodeTreeBuilder {self._tree_type}:{self.name}>"
