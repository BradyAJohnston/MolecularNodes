# Node Builder Design Conclusions

## Core Requirements

1. **IDE-Friendly**: Full autocomplete with type hints for all node types
2. **Scalable**: Must handle hundreds of Blender geometry nodes without manual boilerplate
3. **Readable API**: Clean, intuitive syntax for building node trees
4. **Type-Safe**: Proper type checking for node inputs/outputs

## Initial Implementation (builder.py)

Current approach in `molecularnodes/nodes/builder.py`:
- Manual class definitions for each node type (SetPosition, TransformGeometry, etc.)
- NodeAdder class with individual methods for each node
- Method chaining via node-specific methods (`.set_position()`, `.transform_geometry()`)

### Problems at Scale
- Manual class creation doesn't scale to 100+ nodes
- NodeAdder would need hundreds of methods
- Can't add hundreds of convenience methods to NodeBuilder base class
- Dynamic dispatch (`__getattr__`) breaks IDE autocomplete

## Solution: Code Generation

**Key Insight**: To get IDE autocomplete, you need actual Python classes in actual .py files that IDEs can parse.

### Approach
1. Write a generator script that introspects Blender's node registry
2. Generate Python files with proper class definitions, type hints, and docstrings
3. Check generated code into version control
4. IDEs see everything, autocomplete works perfectly
5. Regenerate when Blender adds new nodes

### Benefits
- ✅ Full IDE autocomplete - real classes in real files
- ✅ Type checking with mypy/pyright
- ✅ Docstrings show in IDE tooltips
- ✅ Maintainable - just regenerate when needed
- ✅ No runtime magic - simple, debuggable code

## Proposed API Design

### Using >> Operator for Linking

```python
from molecularnodes.nodes import tree, Position, SetPosition, TransformGeometry, Vector

with tree.new("MyTree") as t:
    pos = Position()
    offset = Vector(value=(1, 2, 3))

    (
        t.input
        >> SetPosition(position=pos, offset=offset)
        >> TransformGeometry(translation=(0, 0, 1))
        >> t.output
    )
```

### Generated Class Structure

Each node gets a generated class like:

```python
class SetPosition(NodeBuilder):
    """Set the position of points in the geometry."""
    bl_idname = "GeometryNodeSetPosition"

    def __init__(
        self,
        *,
        selection: Boolean | bool | None = None,
        position: Vector3 | tuple[float, float, float] | None = None,
        offset: Vector3 | tuple[float, float, float] | None = None,
    ):
        super().__init__()
        self._establish_links(
            selection=selection,
            position=position,
            offset=offset
        )

    @property
    def geometry(self) -> GeometrySocket:
        """Output geometry"""
        return GeometrySocket(self.node.outputs["Geometry"])
```

### >> Operator Implementation

```python
class NodeBuilder:
    def __rshift__(self, other: "NodeBuilder") -> "NodeBuilder":
        """Chain nodes: self >> other links geometry output to geometry input"""
        # Smart linking: prefer Geometry sockets, fall back to defaults
        self_out = (self.node.outputs.get("Geometry") or
                   self.default_output)
        other_in = (other.node.inputs.get("Geometry") or
                   other.default_input)

        self.tree.link(self_out, other_in)
        return other  # Return target for continued chaining
```

## File Organization

```
molecularnodes/nodes/
├── builder.py          # Core NodeBuilder, TreeBuilder classes (existing)
├── generator.py        # Code generation script (NEW - run manually)
├── _generated/         # Generated code (NEW - checked into git)
│   ├── __init__.py    # Exports all nodes
│   ├── geometry.py    # SetPosition, TransformGeometry, etc.
│   ├── mesh.py        # Mesh nodes
│   ├── curve.py       # Curve nodes
│   ├── attribute.py   # Attribute nodes
│   └── utilities.py   # Math, conversion nodes
└── sockets.py         # Typed socket wrapper classes (NEW)
```

## Additional Features to Add

### Context Manager for Cleaner API
```python
class TreeBuilder:
    _active_tree: ClassVar[TreeBuilder | None] = None

    def __enter__(self):
        TreeBuilder._active_tree = self
        return self

    def __exit__(self, *args):
        TreeBuilder._active_tree = None

class NodeBuilder:
    def __init__(self):
        self._tree = TreeBuilder._active_tree
        if self._tree is None:
            raise RuntimeError("NodeBuilder must be used within tree context")
        self.node = self._tree.add(self.bl_idname)
```

This eliminates needing to pass `tree` to every node.

### Socket Type Wrappers
```python
class GeometrySocket:
    """Strongly typed geometry socket"""
    def __init__(self, socket: NodeSocket):
        self.socket = socket

    def __rshift__(self, other: "NodeBuilder | GeometrySocket"):
        # Can pipe socket directly to another node
        ...
```

Enables type checking: "You can't connect a Boolean to a Vector input"

### Socket Caching
- Build normalized socket name → socket index mapping at node creation
- Cache lookups to avoid repeated string manipulation (lines 157-159, 170-172 in current code)

## Implementation Plan

1. Keep existing NodeBuilder foundation (it's solid)
2. Add >> operator overloading to NodeBuilder
3. Write code generator that introspects Blender nodes
4. Generate initial set of nodes (geometry, mesh, utilities)
5. Add socket wrapper classes for type safety
6. Add context manager for cleaner tree building
7. Add socket caching for performance
8. Document regeneration workflow

## What to Keep from Current Implementation

- TreeBuilder class structure ✅
- NodeBuilder base class ✅
- `_establish_links()` method ✅
- Flexible input handling (values vs connections) ✅
- `link_from()` / `link_to()` methods ✅

## What Was Replaced

- ✅ Manual node class definitions → Will generate from Blender registry
- ✅ NodeAdder class → Removed, use direct imports with context manager
- ✅ Method chaining via `.set_position()` → Now use >> operator
- String type hints → Will use proper TypeAlias definitions (pending)

## Implementation Status

### Completed ✅
1. ✅ Implement >> operator on NodeBuilder (builder.py:~245-263)
2. ✅ Add context manager for automatic tree tracking (builder.py:154-159)
3. ✅ Update NodeBuilder to use active tree from context (builder.py:~196-218)
4. ✅ **Implement named socket access** (builder.py:75-134)
   - Created SocketAccessor class for tree.inputs.socket_name and tree.outputs.socket_name
   - Created SocketNodeBuilder for specific socket access
   - Added name normalization: "Geometry" -> geometry, "My Socket" -> my_socket
   - Plural names (inputs/outputs) match interface(inputs=..., outputs=...) parameters
5. ✅ Update node classes to work with context manager (SetPosition, TransformGeometry, Position, Vector)
6. ✅ Remove old API completely (NodeAdder, method chaining, etc.)
7. ✅ **Implement typed socket interface API** (builder.py:205-255, sockets.py)
   - Created socket definition classes with Socket prefix (SocketGeometry, SocketBoolean, SocketVector, SocketFloat, etc.)
   - Added TreeBuilder.interface() method with inputs/outputs parameters
   - Full IDE autocomplete and type checking for socket definitions
   - Socket prefix disambiguates from node classes (SocketVector vs Vector node)
8. ✅ Create clean examples demonstrating the new API (builder.py:~440-504)

### Next Tasks
1. Prototype the code generator (`generator.py`)
2. Create socket wrapper classes for type safety
3. Generate first batch of node classes from Blender registry
4. Test the new API with existing MolecularNodes use cases
5. Add socket caching for performance optimization
