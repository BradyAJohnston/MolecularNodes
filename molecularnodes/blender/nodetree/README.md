# Node Tree Builder

A fluent, builder-pattern API for simplified Blender node tree creation.

## Overview

This module reduces the boilerplate required when creating Blender node trees by providing:

- **Fluent/chainable API** - Method chaining for concise code
- **Type hints** - Better IDE autocomplete and type checking
- **Smart defaults** - Less configuration for common cases
- **Readable syntax** - Clear intent, less noise

## Quick Start

```python
from molecularnodes.blender.nodetree import NodeTree

# Create a geometry node tree
tree = NodeTree.geometry("My Nodes")

# Add nodes
separate = tree.add_node("ShaderNodeSeparateXYZ", location=(0, 0))
compare = tree.add_node("FunctionNodeCompare",
    name="Compare X",
    operation="LESS_EQUAL",
    location=(200, 0)
)

# Link nodes
separate.output("X").link_to(compare.input(0))

# Auto-arrange and build
tree.auto_layout()
bpy_tree = tree.build()
```

## Benefits

Compared to raw `bpy` API:

**Before (raw bpy):**
```python
tree = bpy.data.node_groups.new("My Nodes", "GeometryNodeTree")
input_node = tree.nodes.new("NodeGroupInput")
output_node = tree.nodes.new("NodeGroupOutput")
input_node.location.x = -200
output_node.location.x = 200
tree.interface.new_socket("Geometry", in_out="INPUT", socket_type="NodeSocketGeometry")
tree.interface.new_socket("Geometry", in_out="OUTPUT", socket_type="NodeSocketGeometry")
tree.links.new(output_node.inputs[0], input_node.outputs[0])

compare = tree.nodes.new("FunctionNodeCompare")
compare.name = "Compare X"
compare.operation = "LESS_EQUAL"
compare.location = (100, 200)
```

**After (fluent API):**
```python
tree = NodeTree.geometry("My Nodes")

compare = tree.add_node("FunctionNodeCompare",
    name="Compare X",
    operation="LESS_EQUAL",
    location=(100, 200)
)
```

**~50% fewer lines, clearer intent**

## Key Features

### 1. Tree Creation

```python
# Geometry node tree (most common)
tree = NodeTree.geometry("My Nodes")

# Shader node tree
tree = NodeTree.shader("My Shader")

# Compositor node tree
tree = NodeTree.compositor("My Compositor")

# Generic tree with manual setup
tree = NodeTree.create("My Tree").init_tree()
```

### 2. Node Creation

```python
# Basic node
node = tree.add_node("FunctionNodeCompare")

# With properties
node = tree.add_node("FunctionNodeCompare",
    name="Compare X",
    operation="LESS_EQUAL",
    location=(100, 200)
)

# Fluent configuration
node = (tree.add_node("ShaderNodeMath")
    .set_name("Multiply")
    .set_property("operation", "MULTIPLY")
    .at(200, 100)
    .set_input_value(1, -1.0)
)

# Custom node groups (GeometryNodeGroup)
style = tree.add_group("Style Density Surface",
    location=(400, 0),
    material="default"
)
```

### 3. Socket Management

```python
# Add individual inputs
tree.add_input("Threshold", "FLOAT",
    description="Comparison threshold",
    default_value=0.5,
    min_value=0.0,
    max_value=1.0
)

# Add multiple inputs at once
tree.add_inputs({
    "Volume": ("GEOMETRY", "Input volume"),
    "ISO Value": ("FLOAT", {"default": 0.8, "min": 0, "max": 1}),
    "Visible": ("BOOL", True),
    "Color": ("COLOR", (1, 0, 0, 1)),
})

# Access sockets
tree.inputs.Geometry  # Attribute access
tree.inputs["Geometry"]  # Dict access
tree.inputs[0]  # Index access
```

### 4. Linking Nodes

```python
# Method 1: Fluent socket API
separate.output("X").link_to(compare.input(0))

# Method 2: Operator syntax
compare.input(0) << separate.output("X")

# Method 3: Tree-level method
tree.link(separate.outputs.X, compare.inputs[0])

# Connect a chain
tree.connect_chain([input_node, process1, process2, output_node])

# Chained linking
(separate.output("X")
    .link_to(compare1.input(0))
    .link_to(compare2.input(0))
)
```

### 5. Node Positioning

```python
# During creation
node = tree.add_node("GeometryNodeJoinGeometry", location=(100, 200))

# Fluent positioning
node.at(100, 200)

# Relative offset
node.offset(50, 0)

# Direct assignment
node.location = (100, 200)

# Auto-layout
tree.auto_layout(spacing=50)
```

### 6. Socket Access

```python
# Access node sockets
node.input("Geometry")  # By name
node.input(0)  # By index
node.inputs.Geometry  # Attribute access
node.inputs["Geometry"]  # Dict access

# Same for outputs
node.output("Result")
node.outputs[0]
node.outputs.Result
```

## Examples

See the `examples/` directory for comprehensive examples:

- `basic_usage.py` - Fundamental API features
- `migration_examples.py` - Before/after comparisons with existing MN code

## Architecture

```
molecularnodes/blender/nodetree/
├── __init__.py          # Public API
├── tree.py              # NodeTreeBuilder class
├── node.py              # NodeWrapper classes
├── socket.py            # SocketWrapper classes
├── DESIGN.md            # Detailed design document
├── README.md            # This file
└── examples/            # Usage examples
    ├── basic_usage.py
    └── migration_examples.py
```

## Type Support

The API is fully typed for IDE autocomplete:

```python
tree: NodeTreeBuilder = NodeTree.geometry("My Nodes")
node: NodeWrapper = tree.add_node("FunctionNodeCompare")
socket: SocketWrapper = node.output("Result")
bpy_tree: bpy.types.NodeTree = tree.build()
```

## Compatibility

- **Fully compatible** with existing code - `tree.build()` returns `bpy.types.NodeTree`
- **Gradual migration** - Mix old and new approaches
- **No performance overhead** - Thin wrapper around bpy API

## API Reference

### `NodeTree` (NodeTreeBuilder)

Main entry point for creating node trees.

**Class Methods:**
- `NodeTree.geometry(name, fallback, input_name, output_name)` - Create geometry tree
- `NodeTree.shader(name)` - Create shader tree
- `NodeTree.compositor(name)` - Create compositor tree
- `NodeTree.create(name)` - Generic tree builder

**Instance Methods:**
- `add_node(type, name, location, **props)` - Add a node
- `add_group(name, location, material, ...)` - Add custom node group
- `add_input(name, type, description, default, **kwargs)` - Add input socket
- `add_output(name, type, description, **kwargs)` - Add output socket
- `add_inputs(dict)` - Add multiple inputs
- `link(from_socket, to_socket)` - Link two sockets
- `connect_chain(nodes)` - Connect nodes in sequence
- `auto_layout(spacing)` - Auto-arrange nodes
- `get_node(name)` - Get node by name
- `remove_node(node)` - Remove a node
- `build()` - Return bpy.types.NodeTree

**Properties:**
- `tree` - The underlying bpy.types.NodeTree
- `input_node` - Group Input node (if exists)
- `output_node` - Group Output node (if exists)
- `inputs` - Input sockets collection
- `outputs` - Output sockets collection

### `NodeWrapper`

Wrapper around bpy.types.Node with fluent API.

**Methods:**
- `set_name(name)` - Set node name
- `set_label(label)` - Set node label
- `set_location(x, y)` / `at(x, y)` - Set position
- `offset(x, y)` - Relative positioning
- `set_property(name, value)` - Set node property
- `set_properties(**kwargs)` - Set multiple properties
- `set_input_value(key, value)` - Set input default value
- `set_input_values(**kwargs)` - Set multiple input values
- `set_parent(parent)` - Set parent node
- `input(key)` - Get input socket
- `output(key)` - Get output socket
- `build()` - Return bpy.types.Node

**Properties:**
- `node` - The underlying bpy.types.Node
- `tree` - The node tree
- `inputs` - Input sockets collection
- `outputs` - Output sockets collection
- `name`, `label`, `location`, `width`

### `GeometryNodeGroupWrapper`

Specialized wrapper for GeometryNodeGroup (extends NodeWrapper).

**Additional Methods:**
- `set_node_tree(tree)` - Assign node tree
- `copy_tree()` - Make tree single-user
- `set_subtree_name(name)` - Name the node tree
- `with_material(material)` - Assign material

### `SocketWrapper`

Wrapper around bpy.types.NodeSocket with linking helpers.

**Methods:**
- `link_to(target)` - Link to target socket
- `unlink()` - Remove all links
- `__lshift__(other)` - Operator: `input << output`

**Properties:**
- `socket` - The underlying bpy.types.NodeSocket
- `node` - The owning node
- `tree` - The node tree
- `name`, `identifier`, `is_output`, `is_linked`, `default_value`

## Development Status

**Current Version:** 0.1.0 (Phase 1 - Fluent Builder API)

**Implemented:**
- ✅ Core fluent API for trees, nodes, sockets
- ✅ Type hints and documentation
- ✅ Examples and migration guide

**Planned:**
- ⏳ Phase 2: Reusable pattern classes (ColorCommon, AnimateFrames, etc.)
- ⏳ Phase 3: Config-based tree generation (dict/YAML)
- ⏳ Unit tests
- ⏳ Integration with existing MN codebase

## Contributing

See `DESIGN.md` for detailed architecture and design decisions.

## License

Part of MolecularNodes - same license as parent project.
