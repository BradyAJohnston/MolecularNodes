# Node Tree Builder - Design Document

## Overview

This module provides a simplified, fluent API for creating Blender node trees, specifically designed to reduce boilerplate and improve readability when working with Geometry Nodes, Shader Nodes, and Compositor Nodes.

## Problem Statement

Creating node trees with the raw `bpy` API is verbose and repetitive:

```python
# Current approach - very verbose (15+ lines for simple setup)
tree = bpy.data.node_groups.new(name="My Nodes", type="GeometryNodeTree")
input_node = tree.nodes.new("NodeGroupInput")
output_node = tree.nodes.new("NodeGroupOutput")
input_node.location.x = -200
output_node.location.x = 200
tree.interface.new_socket("Geometry", in_out="INPUT", socket_type="NodeSocketGeometry")
tree.interface.new_socket("Geometry", in_out="OUTPUT", socket_type="NodeSocketGeometry")
tree.links.new(output_node.inputs[0], input_node.outputs[0])

# Creating and configuring nodes
compare = tree.nodes.new("FunctionNodeCompare")
compare.name = "Compare X Positive"
compare.operation = "LESS_EQUAL"
compare.location = (100, 200)

# Linking
tree.links.new(some_node.outputs["X"], compare.inputs[0])
```

Analysis of MolecularNodes codebase reveals:
- `style_density_iso_surface.py`: ~600 lines for single node group (40+ nodes, 70+ links)
- Repetitive patterns across multiple files
- Manual position management required
- Error-prone socket indexing

## Solution: Hybrid Approach

We implement a **progressive disclosure** API that combines multiple approaches:

### Phase 1: Fluent/Builder Pattern (Foundation)
- Chainable methods for common operations
- Type-safe wrappers around bpy API
- Minimal learning curve
- ~30% reduction in line count

### Phase 2: Reusable Component Classes
- Pre-built patterns common in MolecularNodes:
  - `ColorCommon` pattern
  - `AnimateFrames` pattern
  - Axis comparison/slicing
  - Material assignment
- Encapsulated complex node groups

### Phase 3 (Optional): Declarative Configuration
- Dictionary/YAML based definitions
- Potential for auto-generation from .blend files
- Useful for very large procedural trees

## Architecture

```
molecularnodes/blender/nodetree/
├── DESIGN.md           # This file
├── __init__.py         # Public API exports
├── tree.py             # NodeTree builder class
├── node.py             # Node builder class
├── socket.py           # Socket wrapper with linking helpers
├── patterns/           # Reusable node group patterns (Phase 2)
│   ├── __init__.py
│   ├── color.py
│   ├── animation.py
│   └── comparison.py
└── examples/           # Usage examples
    ├── basic.py
    ├── migration.py
    └── patterns.py
```

## API Design - Phase 1 (Fluent Builder)

### Tree Creation

```python
from molecularnodes.blender.nodetree import NodeTree

# Fluent tree creation
tree = (NodeTree.create("My Geometry Nodes")
    .geometry_tree()
    .add_input("Geometry", "GEOMETRY")
    .add_output("Geometry", "GEOMETRY")
)

# Or more concise for common case:
tree = NodeTree.geometry("My Geometry Nodes")
```

### Node Creation

```python
from molecularnodes.blender.nodetree import Node

# Create node with builder
compare = (Node(tree, "FunctionNodeCompare")
    .set_name("Compare X Positive")
    .set_property("operation", "LESS_EQUAL")
    .set_location(100, 200)
    .build()
)

# Or factory method on tree:
compare = tree.add_node("FunctionNodeCompare",
    name="Compare X Positive",
    operation="LESS_EQUAL",
    location=(100, 200)
)

# For GeometryNodeGroup (custom nodes):
style_node = tree.add_group("Style Density Surface",
    location=(400, 0),
    material="default"
)
```

### Linking

```python
# Method 1: Tree-level linking
tree.link(separate_xyz.outputs["X"], compare.inputs[0])

# Method 2: Fluent socket linking
separate_xyz.output("X").link_to(compare.input(0))

# Method 3: Operator overloading (sugar)
compare.input(0) << separate_xyz.output("X")

# Method 4: Chain connections
(input_node.output(0)
    .link_to(style_node.input(0))
    .link_to(join_node.input(0))
    .link_to(output_node.input(0))
)
```

### Socket Management

```python
# Add multiple inputs at once
tree.add_inputs({
    "volume": ("GEOMETRY", "Input geometry"),
    "iso_value": ("FLOAT", {"min": 0, "max": 1, "default": 0.5}),
    "visible": ("BOOL", True),
    "color": ("COLOR", (0, 0, 1, 1)),
})

# Access sockets via attribute or key
tree.inputs.volume  # Returns socket wrapper
tree.inputs["volume"]  # Alternative syntax
```

### Layout

```python
# Manual positioning
node.at(100, 200)  # Set location
node.offset(-50, 0)  # Relative offset

# Auto-layout (uses existing arrange.py)
tree.auto_layout(spacing=50)
```

## Example Migration

### Before (Current Code)

From `molecularnodes/nodes/nodes.py:300-318`:

```python
def create_starting_nodes_starfile(object):
    node_mod = get_mod(object)
    node_name = f"MN_starfile_{object.name}"

    group = new_tree(node_name)
    node_mod.node_group = group
    link = group.links.new

    node_input = get_input(group)
    node_output = get_output(group)
    node_input.location = [0, 0]
    node_output.location = [700, 0]
    node_star_instances = add_custom(group, "Starfile Instances", [450, 0])
    link(node_star_instances.outputs[0], node_output.inputs[0])
    link(node_input.outputs[0], node_star_instances.inputs[0])
```

### After (With New API)

```python
from molecularnodes.blender.nodetree import NodeTree

def create_starting_nodes_starfile(object):
    node_mod = get_mod(object)
    node_name = f"MN_starfile_{object.name}"

    tree = NodeTree.geometry(node_name)
    node_mod.node_group = tree.build()

    # Fluent chain: input -> custom group -> output
    (tree.input_node.output(0)
        .link_to(tree.add_group("Starfile Instances", location=(450, 0)).input(0))
        .link_to(tree.output_node.input(0))
    )
```

**Result: 50% line reduction, clearer intent**

### More Complex Example

Before (from `nodes.py:321-395`, simplified):

```python
def create_starting_nodes_density(object, threshold=0.8, style="density_surface"):
    mod = get_mod(object)
    node_name = f"MN_density_{object.name}"
    group = new_tree(node_name, fallback=False)
    link = group.links.new
    mod.node_group = group

    node_input = get_input(group)
    node_input.location = [0, 0]
    node_output = get_output(group)
    node_output.location = [800, 0]

    node_density = group.nodes.new("GeometryNodeGroup")
    node_density.name = styles_mapping[style]
    node_density.location = [400, 0]
    tree = style_density_iso_surface_node_group()
    tree.name = f"{styles_mapping[style]}.{object.name}"
    node_density.node_tree = tree
    assign_material(node_density)

    node_join = group.nodes.new("GeometryNodeJoinGeometry")
    node_join.location = [620, 0]

    link(node_input.outputs[0], node_density.inputs[0])
    link(node_density.outputs[0], node_join.inputs[0])
    link(node_join.outputs[0], node_output.inputs[0])

    return node_density
```

After:

```python
def create_starting_nodes_density(object, threshold=0.8, style="density_surface"):
    mod = get_mod(object)
    node_name = f"MN_density_{object.name}"

    tree = NodeTree.geometry(node_name, fallback=False)
    mod.node_group = tree.build()

    # Create and wire nodes in fluent chain
    density = tree.add_group(styles_mapping[style],
        location=(400, 0),
        material="default",
        subtree_name=f"{styles_mapping[style]}.{object.name}"
    )

    join = tree.add_node("GeometryNodeJoinGeometry", location=(620, 0))

    # Wire: input -> density -> join -> output
    tree.connect_chain([tree.input_node, density, join, tree.output_node])

    return density
```

## Design Principles

1. **Progressive Disclosure**: Simple tasks should be simple, complex tasks should be possible
2. **Type Safety**: Use type hints throughout for IDE autocomplete
3. **Minimal Magic**: Stay close to bpy concepts, just reduce boilerplate
4. **Backward Compatible**: Works alongside existing code, allows gradual migration
5. **Fail Fast**: Clear error messages when something goes wrong
6. **Testable**: Pure functions where possible, easy to unit test

## Implementation Priorities

### Phase 1.1 - Core API (Week 1)
- [x] `NodeTree` class with fluent interface
- [ ] `Node` class with fluent interface
- [ ] `Socket` class with linking helpers
- [ ] Basic tree creation (geometry, shader, compositor)
- [ ] Node creation and property setting
- [ ] Link creation
- [ ] Auto-layout integration

### Phase 1.2 - Convenience (Week 2)
- [ ] Batch operations (add_nodes, add_inputs)
- [ ] Smart defaults for common node types
- [ ] Integration with existing `molecularnodes.nodes` helpers
- [ ] Error handling and validation

### Phase 2 - Patterns (Week 3-4)
- [ ] `ColorCommon` pattern
- [ ] `AnimateFrames` pattern
- [ ] Axis comparison patterns
- [ ] Material assignment helpers
- [ ] Common MolecularNodes sub-graphs

### Phase 3 - Advanced (Future)
- [ ] Config-based tree generation
- [ ] Export existing trees to config
- [ ] Visual node graph generator

## Testing Strategy

1. **Unit Tests**: Test each builder method independently
2. **Integration Tests**: Build actual node trees, verify structure
3. **Migration Tests**: Convert existing MN code, compare output
4. **Performance**: Ensure no significant overhead vs raw bpy

## Success Metrics

- [ ] 30-50% reduction in lines of code for typical node tree creation
- [ ] Type hints enable 90%+ IDE autocomplete accuracy
- [ ] Migration guide allows conversion of existing code
- [ ] Zero performance regression
- [ ] Positive developer feedback

## Future Possibilities

1. **Visual Editor Integration**: Generate Python code from Blender UI
2. **Template Library**: Pre-built node graphs for common tasks
3. **Documentation Generation**: Auto-generate docs from node trees
4. **Serialization**: Save/load trees as JSON/YAML
5. **Version Migration**: Handle Blender API changes across versions

## Related Work

- **Blender's bpy**: Low-level API we're wrapping
- **MolecularNodes nodes.py**: Existing helpers we're extending
- **NodeArrange**: Auto-layout system we integrate with
- **Sverchok**: Alternative node creation approach (more visual)

## Questions & Decisions

1. **Q: Should we use operator overloading for links?**
   A: Optional - provide both `link_to()` and `<<` operator, let users choose

2. **Q: How to handle node tree versioning?**
   A: Store builder API version in custom properties for future migration

3. **Q: Integration with existing MN code?**
   A: Make tree.build() return bpy.types.GeometryNodeTree for compatibility

4. **Q: Type hints for socket types?**
   A: Use string literals initially, consider Enums later for stricter typing

## References

- Blender bpy.types.GeometryNodeTree: https://docs.blender.org/api/current/bpy.types.GeometryNodeTree.html
- MolecularNodes codebase analysis: See parent directory analysis
- Builder Pattern: https://refactoring.guru/design-patterns/builder
- Fluent Interface: https://en.wikipedia.org/wiki/Fluent_interface
