# Node Class Generator

This directory contains code to automatically generate Python classes for all Blender geometry nodes.

## Overview

The generator introspects Blender's node registry and creates Python classes with:
- Full type hints for all inputs and outputs
- Automatic parameter handling via `_establish_links()`
- Property accessors for output sockets
- Organized by category (Geometry, Mesh, Curve, etc.)

## How to Run

### Method 1: From Blender GUI

1. Open Blender
2. Go to **Scripting** workspace
3. Open `run_generator.py`
4. Click **Run Script**
5. Check the console for progress
6. Generated files appear in `_generated/` directory

### Method 2: From Command Line

```bash
blender --background --python molecularnodes/nodes/run_generator.py
```

### Method 3: From Blender Python Console

```python
import molecularnodes.nodes.generator as gen
gen.generate_all()
```

## What Gets Generated

### File Structure

```
molecularnodes/nodes/_generated/
├── __init__.py           # Exports all nodes
├── geometry.py           # General geometry nodes
├── mesh.py               # Mesh-specific nodes
├── curve.py              # Curve nodes
├── attribute.py          # Attribute operations
├── input.py              # Input nodes
└── utilities.py          # Utility/math nodes
```

### Example Generated Class

```python
class SetPosition(NodeBuilder):
    """Set the position of every point in the geometry."""
    name = "GeometryNodeSetPosition"

    def __init__(
        self,
        tree: "TreeBuilder | None" = None,
        selection: TYPE_INPUT_BOOLEAN = None,
        position: TYPE_INPUT_VECTOR = None,
        offset: TYPE_INPUT_VECTOR = None
    ):
        super().__init__(tree)
        self._establish_links(
            selection=selection,
            position=position,
            offset=offset
        )

    @property
    def geometry(self) -> NodeSocket:
        """Output socket: Geometry"""
        return self.node.outputs["Geometry"]
```

## How It Works

### 1. Discovery
Iterates through `bpy.types` to find all registered `GeometryNode*` types.

### 2. Introspection
For each node type:
- Creates temporary instance in a temp node tree
- Extracts input sockets (name, type, default values)
- Extracts output sockets (name, type)
- Extracts properties (enums, booleans, etc.)
- Extracts documentation strings

### 3. Code Generation
Generates Python class with:
- `name` attribute (Blender node identifier)
- `__init__` method with all inputs as parameters
- `_establish_links()` call to connect inputs
- `@property` methods for output socket access

### 4. Organization
Groups nodes by category and writes to separate files.

## Inspiration

Based on the approach used by [geometry-script](https://github.com/carson-katri/geometry-script), which dynamically registers nodes at runtime. Our approach differs by:
- Generating static Python files (better IDE support)
- Checking generated code into version control
- Regenerating only when Blender nodes change
- Optimized for type hints and autocomplete

## When to Regenerate

Run the generator when:
- Blender releases a new version with new nodes
- You want to add custom node types
- The generator logic is improved

## Generated Code Quality

**Pros:**
✅ Full IDE autocomplete
✅ Type hints for all parameters
✅ Docstrings from Blender
✅ Property access to outputs
✅ No runtime magic
✅ Easy to debug

**Limitations:**
- Properties and enums need manual refinement for best type hints
- Some complex nodes may need manual adjustments
- Category classification is heuristic-based

## Customization

To customize generation:

1. **Add custom type mappings** in `get_socket_type_hint()`
2. **Improve categorization** in `introspect_node()`
3. **Add custom templates** in `generate_node_class()`
4. **Filter nodes** in `get_all_geometry_nodes()`

## Testing Generated Code

After generation, test with:

```python
from molecularnodes.nodes import TreeBuilder
from molecularnodes.nodes._generated.geometry import SetPosition, TransformGeometry
from molecularnodes.nodes.sockets import SocketGeometry

tree = TreeBuilder("Test")
tree.interface(
    inputs=[SocketGeometry(name="Geometry")],
    outputs=[SocketGeometry(name="Geometry")]
)

with tree:
    tree.inputs.geometry >> SetPosition() >> TransformGeometry() >> tree.outputs.geometry
```

## Troubleshooting

### "Module not found" error
Make sure you're running from within Blender's Python environment.

### No nodes generated
Check Blender version - needs 3.0+ with geometry nodes support.

### Import errors in generated code
The generator may need updates for new Blender API changes.

## Files

- **`generator.py`** - Main generator logic
- **`run_generator.py`** - Convenient runner script
- **`_generated/`** - Output directory (created on first run)

## Future Enhancements

- [ ] Generate enum classes for node properties
- [ ] Add method chaining for common patterns
- [ ] Generate stub files (.pyi) for even better IDE support
- [ ] Support for custom node types
- [ ] Automatic categorization from Blender's menu system
- [ ] Generate tests for each node class
