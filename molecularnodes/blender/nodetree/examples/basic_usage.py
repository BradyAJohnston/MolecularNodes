"""
Basic usage examples for the node tree builder API.

These examples demonstrate the fundamental features of the fluent API.
"""

from molecularnodes.blender.nodetree import NodeTree


def example_simple_geometry_tree():
    """Create a simple geometry node tree with basic nodes."""

    # Create a geometry node tree (includes input/output by default)
    tree = NodeTree.geometry("Simple Example")

    # Add a few nodes with the fluent API
    separate = tree.add_node(
        "ShaderNodeSeparateXYZ",
        name="Separate XYZ",
        location=(0, 0)
    )

    compare = tree.add_node(
        "FunctionNodeCompare",
        name="Compare X",
        operation="LESS_EQUAL",
        location=(200, 0)
    )

    # Link using fluent socket API
    separate.output("X").link_to(compare.input(0))

    # Set input value
    compare.set_input_value(1, 0.5)  # Threshold = 0.5

    # Wire to tree input/output
    tree.inputs.Geometry.link_to(separate.input("Vector"))
    compare.output("Result").link_to(tree.outputs.Geometry)

    # Auto-arrange nodes
    tree.auto_layout()

    # Get the final Blender tree
    return tree.build()


def example_with_custom_sockets():
    """Create a tree with custom input/output sockets."""

    tree = NodeTree.create("Custom Sockets").geometry_tree()

    # Add custom input sockets
    tree.add_input("Threshold", "FLOAT",
        description="Comparison threshold",
        default_value=0.5,
        min_value=0.0,
        max_value=1.0
    )

    tree.add_input("Visible", "BOOL",
        description="Whether geometry is visible",
        default_value=True
    )

    tree.add_input("Color", "COLOR",
        default_value=(1, 0, 0, 1)  # Red
    )

    # Add custom output
    tree.add_output("Selection", "BOOLEAN",
        description="Selection mask"
    )

    return tree.build()


def example_batch_inputs():
    """Add multiple inputs at once using dictionary syntax."""

    tree = NodeTree.geometry("Batch Inputs")

    # Add multiple inputs at once
    tree.add_inputs({
        "Volume": ("GEOMETRY", "Input volume data"),
        "ISO Value": ("FLOAT", {"default": 0.8, "min": 0, "max": 1}),
        "Visible": ("BOOL", True),
        "Positive Color": ("COLOR", (0, 0, 1, 1)),  # Blue
        "Negative Color": ("COLOR", (1, 0, 0, 1)),  # Red
        "Material": "MATERIAL",
    })

    return tree.build()


def example_linking_syntax():
    """Demonstrate different ways to link nodes."""

    tree = NodeTree.geometry("Linking Examples")

    pos = tree.add_node("GeometryNodeInputPosition", location=(0, 0))
    sep = tree.add_node("ShaderNodeSeparateXYZ", location=(200, 0))
    cmp = tree.add_node("FunctionNodeCompare", location=(400, 0))

    # Method 1: Fluent socket linking
    pos.output("Position").link_to(sep.input("Vector"))

    # Method 2: Operator overloading (input << output)
    sep.output("X") >> cmp.input(0)  # Note: using >> for clarity

    # Actually in our API it's:
    cmp.input(0) << sep.output("X")

    # Method 3: Tree-level link method
    tree.link(sep.outputs.Y, cmp.inputs[1])

    return tree.build()


def example_chain_connection():
    """Connect multiple nodes in a chain."""

    tree = NodeTree.geometry("Node Chain")

    # Create several processing nodes
    node1 = tree.add_node("GeometryNodeMeshToPoints", location=(0, 0))
    node2 = tree.add_node("GeometryNodeInstanceOnPoints", location=(200, 0))
    node3 = tree.add_node("GeometryNodeRealizeInstances", location=(400, 0))
    node4 = tree.add_node("GeometryNodeJoinGeometry", location=(600, 0))

    # Connect them all in sequence: input -> n1 -> n2 -> n3 -> n4 -> output
    tree.connect_chain([
        tree.input_node,
        node1,
        node2,
        node3,
        node4,
        tree.output_node
    ])

    return tree.build()


def example_node_properties():
    """Set multiple node properties using different methods."""

    tree = NodeTree.geometry("Node Properties")

    # Method 1: Set properties during creation
    compare = tree.add_node("FunctionNodeCompare",
        name="Compare X Positive",
        operation="LESS_EQUAL",
        data_type="FLOAT",
        location=(100, 200)
    )

    # Method 2: Set properties fluently after creation
    math_node = (tree.add_node("ShaderNodeMath")
        .set_name("Multiply")
        .set_property("operation", "MULTIPLY")
        .set_input_value(1, -1.0)
        .at(200, 100)
    )

    # Method 3: Set multiple properties at once
    boolean = tree.add_node("FunctionNodeBooleanMath")
    boolean.set_properties(
        operation="OR",
        name="Combine XY"
    )

    # Method 4: Set multiple input values
    vector = tree.add_node("FunctionNodeInputVector")
    vector.set_input_values(
        X=1.0,
        Y=2.0,
        Z=3.0
    )

    return tree.build()


def example_node_positioning():
    """Different ways to position nodes."""

    tree = NodeTree.geometry("Positioning")

    # Method 1: Set location during creation
    node1 = tree.add_node("GeometryNodeJoinGeometry", location=(100, 200))

    # Method 2: Set location fluently
    node2 = tree.add_node("GeometryNodeMeshToPoints").at(300, 200)

    # Method 3: Relative offset
    node3 = tree.add_node("GeometryNodeInstanceOnPoints").at(100, 0)
    node3.offset(200, 0)  # Now at (300, 0)

    # Method 4: Set location property directly
    node4 = tree.add_node("GeometryNodeRealizeInstances")
    node4.location = (500, 100)

    return tree.build()


def example_geometry_node_group():
    """Use custom node groups (GeometryNodeGroup)."""

    tree = NodeTree.geometry("Custom Groups")

    # Add a custom MolecularNodes group
    # This loads from the MN data file
    style_node = tree.add_group("Style Density Surface",
        location=(400, 0),
        material="default"
    )

    # Can also set subtree name
    animate_node = tree.add_group("Animate Frames",
        location=(200, 0),
        subtree_name="Animate.MyMolecule"
    )

    # Wire them together
    tree.connect_chain([tree.input_node, animate_node, style_node, tree.output_node])

    return tree.build()


if __name__ == "__main__":
    # Run examples (requires Blender environment)
    print("Running basic usage examples...")

    tree1 = example_simple_geometry_tree()
    print(f"Created: {tree1.name}")

    tree2 = example_with_custom_sockets()
    print(f"Created: {tree2.name}")

    tree3 = example_batch_inputs()
    print(f"Created: {tree3.name}")

    # etc...
