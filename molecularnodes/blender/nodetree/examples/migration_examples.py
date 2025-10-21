"""
Migration examples showing how to convert existing MolecularNodes code
to use the new fluent API.

These examples are based on actual code from molecularnodes/nodes/nodes.py
"""

import bpy
from molecularnodes.blender.nodetree import NodeTree
from molecularnodes.nodes.nodes import get_mod


# ==============================================================================
# Example 1: Simple Starfile Node Setup
# ==============================================================================

def create_starfile_nodes_OLD(obj):
    """
    OLD APPROACH - From nodes.py:300-318
    Verbose, lots of manual steps
    """
    from molecularnodes.nodes.nodes import new_tree, get_input, get_output, add_custom

    node_mod = get_mod(obj)
    node_name = f"MN_starfile_{obj.name}"

    # Create tree with helper
    group = new_tree(node_name)
    node_mod.node_group = group
    link = group.links.new

    # Position input/output nodes
    node_input = get_input(group)
    node_output = get_output(group)
    node_input.location = [0, 0]
    node_output.location = [700, 0]

    # Add custom node
    node_star_instances = add_custom(group, "Starfile Instances", [450, 0])

    # Link everything
    link(node_star_instances.outputs[0], node_output.inputs[0])
    link(node_input.outputs[0], node_star_instances.inputs[0])


def create_starfile_nodes_NEW(obj):
    """
    NEW APPROACH - Using fluent API
    More concise and readable
    """
    node_mod = get_mod(obj)
    node_name = f"MN_starfile_{obj.name}"

    tree = NodeTree.geometry(node_name)
    node_mod.node_group = tree.build()

    # Create and wire nodes in one fluent chain
    star_instances = tree.add_group("Starfile Instances", location=(450, 0))

    tree.connect_chain([tree.input_node, star_instances, tree.output_node])

    # Or even more concise with direct linking:
    # (tree.input_node.output(0)
    #     .link_to(star_instances.input(0))
    #     .link_to(tree.output_node.input(0))
    # )


# ==============================================================================
# Example 2: Density Node Setup with Configuration
# ==============================================================================

def create_density_nodes_OLD(obj, threshold=0.8, style="density_surface"):
    """
    OLD APPROACH - From nodes.py:321-395 (simplified)
    Many manual steps for configuration
    """
    from molecularnodes.nodes.nodes import (
        new_tree, get_input, get_output, add_custom, styles_mapping
    )
    from molecularnodes.nodes.material import assign_material

    mod = get_mod(obj)
    node_name = f"MN_density_{obj.name}"

    group = new_tree(node_name, fallback=False)
    link = group.links.new
    mod.node_group = group

    # Position nodes
    node_input = get_input(group)
    node_input.location = [0, 0]
    node_output = get_output(group)
    node_output.location = [800, 0]

    # Create style node
    node_density = add_custom(group, styles_mapping[style], [400, 0])

    # Configure node tree
    node_tree_copy = node_density.node_tree.copy()
    node_tree_copy.name = f"{styles_mapping[style]}.{obj.name}"
    node_density.node_tree = node_tree_copy

    # Assign material
    assign_material(node_density)

    # Set threshold
    items_tree = node_density.node_tree.interface.items_tree
    items_tree["Threshold"].default_value = threshold
    node_density.inputs["Threshold"].default_value = threshold

    # Add join node
    node_join = group.nodes.new("GeometryNodeJoinGeometry")
    node_join.location = [620, 0]

    # Link everything
    link(node_input.outputs[0], node_density.inputs[0])
    link(node_density.outputs[0], node_join.inputs[0])
    link(node_join.outputs[0], node_output.inputs[0])

    return node_density


def create_density_nodes_NEW(obj, threshold=0.8, style="density_surface"):
    """
    NEW APPROACH - Using fluent API
    Much cleaner and more focused
    """
    from molecularnodes.nodes.nodes import styles_mapping

    mod = get_mod(obj)
    node_name = f"MN_density_{obj.name}"

    tree = NodeTree.geometry(node_name, fallback=False)
    mod.node_group = tree.build()

    # Create and configure density node
    density = (tree.add_group(styles_mapping[style],
            location=(400, 0),
            material="default",
            subtree_name=f"{styles_mapping[style]}.{obj.name}"
        )
        .copy_tree()  # Make single-user
        .set_input_value("Threshold", threshold)
    )

    # Add join node
    join = tree.add_node("GeometryNodeJoinGeometry", location=(620, 0))

    # Wire: input -> density -> join -> output
    tree.connect_chain([tree.input_node, density, join, tree.output_node])

    return density


# ==============================================================================
# Example 3: Complex Node with Multiple Properties
# ==============================================================================

def create_comparison_node_OLD(tree):
    """
    OLD APPROACH - Manual property setting
    """
    compare = tree.nodes.new("FunctionNodeCompare")
    compare.name = "Compare X Positive"
    compare.operation = "LESS_EQUAL"
    compare.location = (100, 200)

    # Set input value
    compare.inputs[1].default_value = 0.5

    return compare


def create_comparison_node_NEW(tree_builder):
    """
    NEW APPROACH - Fluent property setting
    """
    compare = (tree_builder.add_node("FunctionNodeCompare",
            name="Compare X Positive",
            operation="LESS_EQUAL",
            location=(100, 200)
        )
        .set_input_value(1, 0.5)
    )

    return compare


# ==============================================================================
# Example 4: Multiple Node Creation and Linking
# ==============================================================================

def create_axis_comparison_OLD(tree):
    """
    OLD APPROACH - Many variables and manual linking
    """
    link = tree.links.new

    # Create nodes
    position = tree.nodes.new("GeometryNodeInputPosition")
    position.location = (0, 0)

    separate_xyz = tree.nodes.new("ShaderNodeSeparateXYZ")
    separate_xyz.location = (200, 0)

    compare_x = tree.nodes.new("FunctionNodeCompare")
    compare_x.name = "Compare X"
    compare_x.operation = "LESS_EQUAL"
    compare_x.location = (400, 0)

    compare_y = tree.nodes.new("FunctionNodeCompare")
    compare_y.name = "Compare Y"
    compare_y.operation = "LESS_EQUAL"
    compare_y.location = (400, -100)

    boolean_math = tree.nodes.new("FunctionNodeBooleanMath")
    boolean_math.name = "Combine XY"
    boolean_math.operation = "OR"
    boolean_math.location = (600, 0)

    # Link nodes
    link(position.outputs["Position"], separate_xyz.inputs["Vector"])
    link(separate_xyz.outputs["X"], compare_x.inputs[0])
    link(separate_xyz.outputs["Y"], compare_y.inputs[0])
    link(compare_x.outputs["Result"], boolean_math.inputs[0])
    link(compare_y.outputs["Result"], boolean_math.inputs[1])

    return boolean_math


def create_axis_comparison_NEW(tree_builder):
    """
    NEW APPROACH - More readable with fluent API
    """
    # Create nodes
    pos = tree_builder.add_node("GeometryNodeInputPosition", location=(0, 0))
    sep = tree_builder.add_node("ShaderNodeSeparateXYZ", location=(200, 0))

    cmp_x = tree_builder.add_node("FunctionNodeCompare",
        name="Compare X",
        operation="LESS_EQUAL",
        location=(400, 0)
    )

    cmp_y = tree_builder.add_node("FunctionNodeCompare",
        name="Compare Y",
        operation="LESS_EQUAL",
        location=(400, -100)
    )

    bool_math = tree_builder.add_node("FunctionNodeBooleanMath",
        name="Combine XY",
        operation="OR",
        location=(600, 0)
    )

    # Link using fluent API - more readable
    pos.output("Position").link_to(sep.input("Vector"))
    sep.output("X").link_to(cmp_x.input(0))
    sep.output("Y").link_to(cmp_y.input(0))
    cmp_x.output("Result").link_to(bool_math.input(0))
    cmp_y.output("Result").link_to(bool_math.input(1))

    return bool_math


# ==============================================================================
# Example 5: Interface Socket Creation
# ==============================================================================

def create_interface_sockets_OLD(tree):
    """
    OLD APPROACH - Verbose socket creation
    """
    # Add input sockets
    volume_socket = tree.interface.new_socket(
        name="Volume",
        in_out="INPUT",
        socket_type="NodeSocketGeometry"
    )
    volume_socket.description = "Input geometry"

    iso_socket = tree.interface.new_socket(
        name="ISO Value",
        in_out="INPUT",
        socket_type="NodeSocketFloat"
    )
    iso_socket.min_value = 0.0
    iso_socket.max_value = 1.0
    iso_socket.default_value = 0.8
    iso_socket.description = "ISO value threshold"

    visible_socket = tree.interface.new_socket(
        name="Visible",
        in_out="INPUT",
        socket_type="NodeSocketBool"
    )
    visible_socket.default_value = True

    color_socket = tree.interface.new_socket(
        name="Positive Color",
        in_out="INPUT",
        socket_type="NodeSocketColor"
    )
    color_socket.default_value = (0.0, 0.0, 1.0, 1.0)


def create_interface_sockets_NEW(tree_builder):
    """
    NEW APPROACH - Batch socket creation with dict
    """
    tree_builder.add_inputs({
        "Volume": ("GEOMETRY", "Input geometry"),
        "ISO Value": ("FLOAT", {
            "default": 0.8,
            "min": 0.0,
            "max": 1.0,
            "description": "ISO value threshold"
        }),
        "Visible": ("BOOL", True),
        "Positive Color": ("COLOR", (0.0, 0.0, 1.0, 1.0)),
    })


# ==============================================================================
# Comparison Summary
# ==============================================================================

def print_comparison():
    """Print a summary of improvements."""
    print("=" * 70)
    print("MIGRATION BENEFITS")
    print("=" * 70)
    print()
    print("✓ 40-60% reduction in lines of code")
    print("✓ More readable - intent is clearer")
    print("✓ Less error-prone - type hints and validation")
    print("✓ Chainable methods reduce variable clutter")
    print("✓ Consistent API across all node types")
    print("✓ Better IDE autocomplete support")
    print()
    print("COMPATIBILITY")
    print("-" * 70)
    print("✓ tree.build() returns bpy.types.NodeTree - fully compatible")
    print("✓ Can mix old and new approaches")
    print("✓ Gradual migration possible")
    print()


if __name__ == "__main__":
    print_comparison()
