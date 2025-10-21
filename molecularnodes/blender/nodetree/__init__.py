"""
Fluent API for simplified Blender node tree creation.

This module provides a builder pattern API for creating Blender node trees
with less boilerplate and more readable code.

Example:
    from molecularnodes.blender.nodetree import NodeTree

    # Create a geometry node tree
    tree = NodeTree.geometry("My Nodes")

    # Add nodes
    compare = tree.add_node("FunctionNodeCompare",
        name="Compare X",
        operation="LESS_EQUAL",
        location=(100, 200)
    )

    separate = tree.add_node("ShaderNodeSeparateXYZ", location=(0, 200))

    # Link nodes
    separate.output("X").link_to(compare.input(0))

    # Or use operator syntax
    compare.input(0) << separate.output("X")

    # Get the final tree
    bpy_tree = tree.build()
"""

from .tree import NodeTreeBuilder as NodeTree
from .node import NodeWrapper, GeometryNodeGroupWrapper
from .socket import SocketWrapper, SocketCollection

__all__ = [
    "NodeTree",
    "NodeWrapper",
    "GeometryNodeGroupWrapper",
    "SocketWrapper",
    "SocketCollection",
]

__version__ = "0.1.0"
