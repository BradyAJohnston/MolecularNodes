# Socket Access: Discoverability Options

## The Challenge

How do we make socket access discoverable in the IDE?

```python
tree.interface(
    inputs=[
        Geometry(name="Geometry"),
        Boolean(name="Selection"),
        Vector(name="Offset"),
    ],
    outputs=[Geometry(name="Geometry")]
)

with tree:
    tree.input.???  # ← How does the IDE know what's available?
```

With `__getattr__`, the IDE can't autocomplete socket names without additional help.

---

## Option 1: Plural Names for Clarity

Use `inputs` and `outputs` (plural) to match the interface definition.

```python
tree.interface(
    inputs=[...],
    outputs=[...]
)

with tree:
    tree.inputs.geometry     # Plural matches interface parameter name
    tree.outputs.geometry
```

**Pros:**
- Matches the `interface(inputs=..., outputs=...)` parameter names
- More grammatically correct ("inputs" is a collection)
- Clear that you're accessing from a collection

**Cons:**
- Still doesn't solve IDE autocomplete issue
- Slightly longer to type

---

## Option 2: Generate Typed Accessor Class Per Tree

When you define the interface, generate a typed class specifically for that tree.

```python
tree = TreeBuilder("MyTree")

tree.interface(
    inputs=[
        Geometry(name="Geometry"),
        Boolean(name="Selection"),
        Vector(name="Offset"),
    ],
    outputs=[Geometry(name="Geometry")]
)

# Behind the scenes, this generates:
# class MyTreeInputs:
#     @property
#     def geometry(self) -> GeometrySocket: ...
#     @property
#     def selection(self) -> BoolSocket: ...
#     @property
#     def offset(self) -> VectorSocket: ...
#
# tree.inputs = MyTreeInputs(tree)

with tree:
    tree.inputs.geometry  # ← IDE autocomplete works! (at runtime)
    tree.inputs.selection
```

**Implementation sketch:**
```python
def interface(self, inputs, outputs):
    # Create sockets...

    # Generate typed accessor class
    input_props = {}
    for socket_def in inputs:
        normalized = self._normalize_name(socket_def.name)
        socket_type = self._get_socket_type(socket_def)

        def make_prop(name, sock_type):
            @property
            def prop(self) -> sock_type:
                return sock_type(self.tree.input(), name)
            return prop

        input_props[normalized] = make_prop(socket_def.name, socket_type)

    InputAccessor = type(f'{self.tree.name}Inputs', (), input_props)
    self.inputs = InputAccessor()
```

**Pros:**
- Runtime autocomplete (if using iPython/Jupyter)
- Type hints work at runtime
- Each tree gets a custom accessor

**Cons:**
- IDE autocomplete still doesn't work (IDE doesn't see runtime classes)
- Complex implementation
- No static type checking

---

## Option 3: Generate .pyi Stub File

Generate a `.pyi` type stub file for each tree that IDEs can read.

```python
tree = TreeBuilder("MyTree")
tree.interface(inputs=[...], outputs=[...])
tree.generate_stubs()  # Creates my_tree.pyi
```

**Generated my_tree.pyi:**
```python
from molecularnodes.nodes import TreeBuilder
from molecularnodes.nodes.sockets import GeometrySocket, BoolSocket, VectorSocket

class MyTreeInputs:
    @property
    def geometry(self) -> GeometrySocket: ...
    @property
    def selection(self) -> BoolSocket: ...
    @property
    def offset(self) -> VectorSocket: ...

class MyTreeOutputs:
    @property
    def geometry(self) -> GeometrySocket: ...

class MyTree(TreeBuilder):
    inputs: MyTreeInputs
    outputs: MyTreeOutputs
```

**Pros:**
- Full IDE autocomplete
- Static type checking
- Standard Python approach

**Cons:**
- Requires file generation step
- Stubs need to be updated when tree changes
- More complex workflow

---

## Option 4: Pre-Import Common Socket Accessors

For common tree patterns, provide pre-made accessor classes.

```python
from molecularnodes.nodes import TreeBuilder
from molecularnodes.nodes.accessors import SimpleGeometryIO, AtomIO, SurfaceIO

# Pre-defined accessor class
tree = TreeBuilder("MyTree", accessor=SimpleGeometryIO)
# SimpleGeometryIO defines: inputs.geometry, outputs.geometry

with tree:
    tree.inputs.geometry  # ← IDE knows about this!
```

**Pre-defined in molecularnodes/nodes/accessors.py:**
```python
@dataclass
class SimpleGeometryIO:
    """Single geometry input and output."""

    @dataclass
    class Inputs:
        @property
        def geometry(self) -> GeometrySocket:
            """Geometry input socket."""
            ...

    @dataclass
    class Outputs:
        @property
        def geometry(self) -> GeometrySocket:
            """Geometry output socket."""
            ...
```

**Pros:**
- Full IDE autocomplete for common patterns
- Type hints work perfectly
- Simple to use for standard cases

**Cons:**
- Only works for pre-defined patterns
- Not flexible for custom interfaces
- Defeats the purpose of dynamic interface definition

---

## Option 5: Protocol + Comment Hint (Pragmatic)

Use Protocol for type hinting with an inline comment helper.

```python
tree = TreeBuilder("MyTree")
tree.interface(
    inputs=[
        Geometry(name="Geometry"),
        Boolean(name="Selection"),
        Vector(name="Offset"),
    ],
    outputs=[Geometry(name="Geometry")]
)

with tree:
    # Available inputs: geometry, selection, offset
    # Available outputs: geometry

    tree.inputs.geometry
    tree.inputs.selection
```

**With Protocol for some type safety:**
```python
from typing import Protocol

class HasGeometryInput(Protocol):
    @property
    def inputs(self) -> object: ...  # IDE at least knows .inputs exists

tree: HasGeometryInput = TreeBuilder("MyTree")
```

**Pros:**
- Simple, no magic
- Comments provide documentation
- Works today without complex implementation

**Cons:**
- Comments can get out of sync
- No actual IDE autocomplete
- Manual documentation

---

## Option 6: Expose via Index Access + Helper

Keep dictionary-style access but provide a helper to show available sockets.

```python
tree.interface(
    inputs=[
        Geometry(name="Geometry"),
        Boolean(name="Selection"),
    ],
    outputs=[Geometry(name="Geometry")]
)

# Helper to see what's available
tree.list_inputs()   # Prints: geometry, selection
tree.list_outputs()  # Prints: geometry

with tree:
    # Use dictionary access (no autocomplete, but clear)
    tree.inputs["geometry"]
    tree.inputs["selection"]

    # Or use normalized names (some autocomplete via __getattr__)
    tree.inputs.geometry
    tree.inputs.selection
```

**Pros:**
- Explicit helper for discoverability
- Dictionary access as fallback
- Simple implementation

**Cons:**
- Still no IDE autocomplete
- Manual discovery step

---

## Option 7: Decorator Pattern (Most Discoverable!)

Use a decorator to define the tree interface, generating a proper typed class.

```python
from molecularnodes.nodes import node_tree, NodeTreeBuilder
from molecularnodes.nodes.sockets import Geometry, Boolean, Vector

@node_tree(
    inputs=[
        Geometry(name="Geometry"),
        Boolean(name="Selection", default=True),
        Vector(name="Offset", default=(0, 0, 0)),
    ],
    outputs=[
        Geometry(name="Geometry"),
    ]
)
class MyTree(NodeTreeBuilder):
    """Custom node tree with typed inputs/outputs."""

    def build(self):
        # IDE knows about these!
        self.inputs.geometry
        self.inputs.selection
        self.inputs.offset
        self.outputs.geometry
```

**Implementation:**
```python
def node_tree(inputs, outputs):
    def decorator(cls):
        # Generate typed input/output accessor classes
        InputsClass = _generate_accessor_class("Inputs", inputs)
        OutputsClass = _generate_accessor_class("Outputs", outputs)

        # Add to class
        cls.__annotations__['inputs'] = InputsClass
        cls.__annotations__['outputs'] = OutputsClass

        return cls
    return decorator
```

**Pros:**
- **Full IDE autocomplete and type hints**
- Class-based, familiar pattern
- Self-documenting
- Scalable to code generation

**Cons:**
- Different pattern from current TreeBuilder
- More complex decorator implementation
- Requires class definition (but could still generate!)

---

## Recommendation: Hybrid Approach

**For most users (Option 6 Enhanced):**
```python
tree = TreeBuilder("MyTree")
tree.interface(
    inputs=[
        Geometry(name="Geometry"),
        Boolean(name="Selection"),
        Vector(name="Offset"),
    ],
    outputs=[Geometry(name="Geometry")]
)

# Use plural for clarity
with tree:
    tree.inputs.geometry    # Property access (works, some IDE support)
    tree.inputs.selection
```

**For generated/library code (Option 7):**
```python
# Code-generated trees get proper type hints
@node_tree(inputs=[...], outputs=[...])
class ColorByAttribute(NodeTreeBuilder):
    def build(self):
        self.inputs.geometry  # Full IDE autocomplete
```

**Implementation priorities:**
1. Use `inputs`/`outputs` (plural) for consistency
2. Implement `__getattr__` for dynamic access (works now)
3. Add `list_inputs()` / `list_outputs()` helper methods
4. Later: Generate stub files or typed classes for library trees

This gives us:
- ✅ Works immediately with basic __getattr__
- ✅ Plural names match interface definition
- ✅ Discoverable via helper methods
- ✅ Path to full type hints via code generation
- ✅ Consistent with design principles

What do you think? Should we go with `inputs`/`outputs` and implement the basic accessor first?
