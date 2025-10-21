"""
Node wrapper providing fluent API for Blender node creation and configuration.
"""

from typing import TYPE_CHECKING, Any, Union
from mathutils import Vector

import bpy

from .socket import SocketWrapper, SocketCollection

if TYPE_CHECKING:
    from .tree import NodeTreeBuilder


class NodeWrapper:
    """
    Wrapper around bpy.types.Node providing fluent configuration API.

    Examples:
        # Create and configure node
        node = (NodeWrapper(tree, "FunctionNodeCompare")
            .set_name("Compare X")
            .set_property("operation", "LESS_EQUAL")
            .set_location(100, 200)
            .build()
        )

        # Access sockets
        node.output("Result").link_to(other_node.input(0))
        node.outputs.Result.link_to(other_node.inputs[0])
    """

    def __init__(
        self,
        tree: Union["NodeTreeBuilder", bpy.types.NodeTree],
        node_type: str,
        node: bpy.types.Node | None = None,
    ):
        """
        Initialize node wrapper.

        Args:
            tree: The node tree (can be NodeTreeBuilder or raw bpy tree)
            node_type: Blender node type identifier (e.g., "FunctionNodeCompare")
            node: Existing node to wrap (if None, creates new node)
        """
        # Handle both NodeTreeBuilder and raw bpy.types.NodeTree
        if hasattr(tree, "_tree"):
            self._tree_builder = tree
            self._tree = tree._tree
        else:
            self._tree_builder = None
            self._tree = tree

        if node is None:
            # Create new node
            self._node = self._tree.nodes.new(node_type)
        else:
            # Wrap existing node
            self._node = node

        self._node_type = node_type

    @property
    def node(self) -> bpy.types.Node:
        """Get the underlying Blender node."""
        return self._node

    @property
    def tree(self) -> bpy.types.NodeTree:
        """Get the node tree this node belongs to."""
        return self._tree

    @property
    def name(self) -> str:
        """Get the node name."""
        return self._node.name

    @name.setter
    def name(self, value: str):
        """Set the node name."""
        self._node.name = value

    @property
    def label(self) -> str:
        """Get the node label."""
        return self._node.label

    @label.setter
    def label(self, value: str):
        """Set the node label."""
        self._node.label = value

    @property
    def location(self) -> Vector:
        """Get the node location."""
        return self._node.location

    @location.setter
    def location(self, value: tuple[float, float] | Vector):
        """Set the node location."""
        self._node.location = value

    @property
    def width(self) -> float:
        """Get the node width."""
        return self._node.width

    @width.setter
    def width(self, value: float):
        """Set the node width."""
        self._node.width = value

    @property
    def inputs(self) -> SocketCollection:
        """Get the node's input sockets."""
        return SocketCollection(self._node.inputs)

    @property
    def outputs(self) -> SocketCollection:
        """Get the node's output sockets."""
        return SocketCollection(self._node.outputs)

    def input(self, key: Union[str, int]) -> SocketWrapper:
        """
        Get an input socket by name or index.

        Args:
            key: Socket name or index

        Returns:
            Wrapped input socket
        """
        return SocketWrapper(self._node.inputs[key])

    def output(self, key: Union[str, int]) -> SocketWrapper:
        """
        Get an output socket by name or index.

        Args:
            key: Socket name or index

        Returns:
            Wrapped output socket
        """
        return SocketWrapper(self._node.outputs[key])

    def set_name(self, name: str) -> "NodeWrapper":
        """
        Set the node name (fluent).

        Args:
            name: New node name

        Returns:
            Self for chaining
        """
        self._node.name = name
        return self

    def set_label(self, label: str) -> "NodeWrapper":
        """
        Set the node label (fluent).

        Args:
            label: New node label

        Returns:
            Self for chaining
        """
        self._node.label = label
        return self

    def set_location(self, x: float, y: float) -> "NodeWrapper":
        """
        Set the node location (fluent).

        Args:
            x: X coordinate
            y: Y coordinate

        Returns:
            Self for chaining
        """
        self._node.location = (x, y)
        return self

    def at(self, x: float, y: float) -> "NodeWrapper":
        """
        Alias for set_location (more concise).

        Args:
            x: X coordinate
            y: Y coordinate

        Returns:
            Self for chaining
        """
        return self.set_location(x, y)

    def offset(self, x: float, y: float) -> "NodeWrapper":
        """
        Offset the node location relative to current position.

        Args:
            x: X offset
            y: Y offset

        Returns:
            Self for chaining
        """
        self._node.location += Vector((x, y))
        return self

    def set_width(self, width: float) -> "NodeWrapper":
        """
        Set the node width (fluent).

        Args:
            width: New node width

        Returns:
            Self for chaining
        """
        self._node.width = width
        return self

    def set_property(self, name: str, value: Any) -> "NodeWrapper":
        """
        Set a node property (fluent).

        Args:
            name: Property name (e.g., "operation", "data_type")
            value: Property value

        Returns:
            Self for chaining
        """
        setattr(self._node, name, value)
        return self

    def set_properties(self, **kwargs) -> "NodeWrapper":
        """
        Set multiple node properties at once.

        Args:
            **kwargs: Property name-value pairs

        Returns:
            Self for chaining

        Example:
            node.set_properties(
                operation="LESS_EQUAL",
                data_type="FLOAT"
            )
        """
        for name, value in kwargs.items():
            setattr(self._node, name, value)
        return self

    def set_input_value(self, key: Union[str, int], value: Any) -> "NodeWrapper":
        """
        Set an input socket's default value.

        Args:
            key: Socket name or index
            value: Default value to set

        Returns:
            Self for chaining
        """
        self._node.inputs[key].default_value = value
        return self

    def set_input_values(self, **kwargs) -> "NodeWrapper":
        """
        Set multiple input socket default values.

        Args:
            **kwargs: Socket name-value pairs

        Returns:
            Self for chaining

        Example:
            node.set_input_values(
                Threshold=0.5,
                Scale=2.0
            )
        """
        for name, value in kwargs.items():
            self._node.inputs[name].default_value = value
        return self

    def set_parent(self, parent: Union["NodeWrapper", bpy.types.Node]) -> "NodeWrapper":
        """
        Set the parent node (for frames).

        Args:
            parent: Parent node (wrapper or raw node)

        Returns:
            Self for chaining
        """
        parent_node = parent.node if isinstance(parent, NodeWrapper) else parent
        self._node.parent = parent_node
        return self

    def build(self) -> bpy.types.Node:
        """
        Return the underlying Blender node (finalizes builder).

        Returns:
            The raw bpy.types.Node
        """
        return self._node

    def __repr__(self) -> str:
        return f"<NodeWrapper {self._node_type}:{self.name}>"


class GeometryNodeGroupWrapper(NodeWrapper):
    """
    Specialized wrapper for GeometryNodeGroup nodes (custom node groups).

    Adds convenience methods for node tree assignment and material handling.
    """

    def __init__(
        self,
        tree: Union["NodeTreeBuilder", bpy.types.NodeTree],
        node: bpy.types.GeometryNodeGroup | None = None,
    ):
        """
        Initialize geometry node group wrapper.

        Args:
            tree: The node tree
            node: Existing GeometryNodeGroup node (if None, creates new)
        """
        super().__init__(tree, "GeometryNodeGroup", node)

    @property
    def node_tree(self) -> bpy.types.NodeTree | None:
        """Get the node tree assigned to this group."""
        return self._node.node_tree

    @node_tree.setter
    def node_tree(self, tree: bpy.types.NodeTree):
        """Set the node tree assigned to this group."""
        self._node.node_tree = tree

    def set_node_tree(self, tree: bpy.types.NodeTree) -> "GeometryNodeGroupWrapper":
        """
        Set the node tree (fluent).

        Args:
            tree: Node tree to assign

        Returns:
            Self for chaining
        """
        self._node.node_tree = tree
        return self

    def copy_tree(self) -> "GeometryNodeGroupWrapper":
        """
        Make the node tree single-user (copy it).

        Returns:
            Self for chaining
        """
        if self._node.node_tree:
            self._node.node_tree = self._node.node_tree.copy()
        return self

    def set_subtree_name(self, name: str) -> "GeometryNodeGroupWrapper":
        """
        Set the name of the assigned node tree.

        Args:
            name: New tree name

        Returns:
            Self for chaining
        """
        if self._node.node_tree:
            self._node.node_tree.name = name
        return self

    def with_material(self, material: str | bpy.types.Material) -> "GeometryNodeGroupWrapper":
        """
        Assign a material to the Material input (if it exists).

        Args:
            material: Material name or material object

        Returns:
            Self for chaining
        """
        # Import here to avoid circular dependency
        from ...nodes.material import assign_material

        assign_material(self._node, new_material=material)
        return self

    def __repr__(self) -> str:
        tree_name = self.node_tree.name if self.node_tree else "None"
        return f"<GeometryNodeGroupWrapper {self.name} tree:{tree_name}>"
