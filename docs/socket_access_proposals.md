# Socket Access Proposals

## Current Problem

With the new typed interface:

```python
tree.interface(
    inputs=[
        Geometry(name="Geometry"),
        Boolean(name="Selection", default=True),
        Vector(name="Offset", default=(0, 0, 0)),
    ],
    outputs=[
        Geometry(name="Geometry"),
        Integer(name="Count"),
    ]
)
```

How do we access these specific sockets in the node chain?

**Current approach is limited:**
```python
tree.input_socket >> SetPosition() >> tree.output_socket
# Problem: Which input socket? Which output socket?
```

---

## Proposal 1: Named Property Access (Recommended)

Access sockets by normalized property names.

```python
tree.interface(
    inputs=[
        Geometry(name="Geometry"),
        Boolean(name="Selection", default=True),
        Vector(name="Offset", default=(0, 0, 0)),
    ],
    outputs=[
        Geometry(name="Geometry"),
    ]
)

with tree:
    (
        tree.input.geometry  # Access "Geometry" input by name
        >> SetPosition(
            selection=tree.input.selection,  # Access "Selection" input
            offset=tree.input.offset          # Access "Offset" input
        )
        >> tree.output.geometry  # Access "Geometry" output by name
    )
```

**Implementation:**
```python
class SocketAccessor:
    """Provides named access to interface sockets."""

    def __init__(self, tree: TreeBuilder, direction: str):
        self.tree = tree
        self.direction = direction  # 'INPUT' or 'OUTPUT'

    def __getattr__(self, name: str):
        # Normalize name (geometry, selection, offset, etc.)
        socket_name = self._denormalize_name(name)

        if self.direction == 'INPUT':
            node = self.tree.input()
        else:
            node = self.tree.output()

        # Return a socket wrapper that can be used in >> chains
        return SocketOutput(node, socket_name)

    def _denormalize_name(self, attr_name: str) -> str:
        """Convert 'geometry' -> 'Geometry', 'my_socket' -> 'My Socket'"""
        return attr_name.replace('_', ' ').title()

class TreeBuilder:
    def __init__(self, ...):
        self.input = SocketAccessor(self, 'INPUT')
        self.output = SocketAccessor(self, 'OUTPUT')
```

**Pros:**
- Clean, readable syntax
- IDE autocomplete (with __getattr__ or dynamic properties)
- Explicit about which socket you're accessing
- Natural Python property access

**Cons:**
- Requires name normalization (spaces to underscores)
- __getattr__ doesn't provide perfect IDE autocomplete (but works at runtime)

---

## Proposal 2: Dictionary-Style Access

Use bracket notation for socket access.

```python
with tree:
    (
        tree.input["Geometry"]
        >> SetPosition(
            selection=tree.input["Selection"],
            offset=tree.input["Offset"]
        )
        >> tree.output["Geometry"]
    )
```

**Implementation:**
```python
class SocketAccessor:
    def __getitem__(self, socket_name: str):
        if self.direction == 'INPUT':
            node = self.tree.input()
        else:
            node = self.tree.output()
        return SocketOutput(node, socket_name)
```

**Pros:**
- No name normalization needed
- Explicit string access
- Familiar Python pattern

**Cons:**
- No IDE autocomplete for socket names
- More verbose (quotes required)
- String-based (defeats our type-safety goals)

---

## Proposal 3: Return Socket References from interface()

Have `interface()` return an object with named attributes.

```python
io = tree.interface(
    inputs=[
        Geometry(name="Geometry"),
        Boolean(name="Selection", default=True),
        Vector(name="Offset", default=(0, 0, 0)),
    ],
    outputs=[
        Geometry(name="Geometry"),
    ]
)

with tree:
    (
        io.inputs.geometry
        >> SetPosition(
            selection=io.inputs.selection,
            offset=io.inputs.offset
        )
        >> io.outputs.geometry
    )
```

**Implementation:**
```python
@dataclass
class InterfaceReference:
    inputs: SocketAccessor
    outputs: SocketAccessor

def interface(self, inputs, outputs):
    # Create sockets...

    return InterfaceReference(
        inputs=SocketAccessor(self, 'INPUT'),
        outputs=SocketAccessor(self, 'OUTPUT')
    )
```

**Pros:**
- Explicit interface reference
- Could enable multiple interface definitions
- Clear separation from tree object

**Cons:**
- Extra variable to track
- Less convenient than using `tree.input`/`tree.output`

---

## Proposal 4: Generate Typed Accessor Classes

Use the interface definition to generate a typed accessor at runtime.

```python
tree = TreeBuilder("MyTree")

tree.interface(
    inputs=[
        Geometry(name="Geometry"),
        Boolean(name="Selection", default=True),
    ],
    outputs=[
        Geometry(name="Geometry"),
    ]
)

# Behind the scenes, generates:
# class MyTreeInputs:
#     @property
#     def geometry(self) -> GeometrySocket: ...
#     @property
#     def selection(self) -> BoolSocket: ...

with tree:
    tree.input.geometry      # Full type hints!
    tree.input.selection     # IDE knows this is BoolSocket!
```

**Implementation:**
```python
def interface(self, inputs, outputs):
    # Create sockets...

    # Dynamically generate typed accessor
    input_properties = {}
    for socket_def in inputs:
        normalized_name = self._normalize_name(socket_def.name)
        socket_type = self._get_socket_wrapper_type(socket_def)

        # Create property dynamically
        def make_getter(socket_name):
            def getter(self):
                return SocketOutput(self.tree.input(), socket_name)
            return property(getter)

        input_properties[normalized_name] = make_getter(socket_def.name)

    # Create class dynamically
    InputAccessor = type('InputAccessor', (), input_properties)
    self.input = InputAccessor()
```

**Pros:**
- **Full type hints and IDE autocomplete** (if done right)
- Compile-time type checking
- Best developer experience

**Cons:**
- Complex implementation
- Runtime class generation is advanced Python
- Harder to debug

---

## Proposal 5: Hybrid - Smart Defaults + Named Access

Use `input_socket`/`output_socket` for simple cases, named access for complex ones.

```python
# Simple case - only one geometry in/out
tree.interface(
    inputs=[Geometry(name="Geometry")],
    outputs=[Geometry(name="Geometry")]
)

with tree:
    tree.input_socket >> SetPosition() >> tree.output_socket

# Complex case - named access
tree.interface(
    inputs=[
        Geometry(name="Atoms"),
        Geometry(name="Surface"),
        Boolean(name="Selection"),
    ],
    outputs=[
        Geometry(name="Result"),
    ]
)

with tree:
    tree.input.atoms >> SetPosition() >> tree.output.result
    tree.input.surface >> MergeGeometry() >> tree.output.result
```

**Pros:**
- Progressive disclosure of complexity
- Simple cases stay simple
- Flexible for complex cases

**Cons:**
- Two different patterns to learn
- Inconsistent API

---

## Recommendation: Proposal 1 + Enhancement

Use **named property access** with a helper for better IDE support:

```python
tree.interface(
    inputs=[
        Geometry(name="Geometry"),
        Boolean(name="Selection", default=True),
        Vector(name="Offset", default=(0, 0, 0)),
    ],
    outputs=[
        Geometry(name="Geometry"),
    ]
)

with tree:
    (
        tree.input.geometry
        >> SetPosition(
            selection=tree.input.selection,
            offset=tree.input.offset
        )
        >> tree.output.geometry
    )
```

### Enhanced with Socket Wrapper Classes

Return typed socket wrappers that know their type:

```python
class GeometrySocket:
    """Wrapper for a geometry socket with >> operator support."""
    def __init__(self, node: Node, socket_name: str):
        self.node = node
        self.socket_name = socket_name
        self.socket = node.outputs[socket_name]

    def __rshift__(self, other: NodeBuilder) -> NodeBuilder:
        # Connect this socket to the other node
        tree = TreeBuilder._active_tree
        tree.link(self.socket, other.default_input)
        return other
```

This gives us:
1. ✅ Clean syntax
2. ✅ Named access to sockets
3. ✅ Works with >> operator
4. ✅ Can be typed for better IDE support
5. ✅ Consistent with typed interface design

---

## Additional Consideration: Socket Outputs from Nodes

What about accessing specific outputs from nodes?

```python
# Current - assumes default output
pos = Position()
tree.input.geometry >> SetPosition(position=pos)

# What if Position has multiple outputs?
pos = Position()
tree.input.geometry >> SetPosition(position=pos.position)  # Specific output

# Or for nodes with multiple geometry outputs?
separate = SeparateGeometry()
separate.selection >> other_node
separate.inverted >> another_node
```

This same named-access pattern could apply to node outputs!

```python
class Position(NodeBuilder):
    @property
    def position(self) -> VectorSocket:
        return VectorSocket(self.node, "Position")
```

Already partially implemented in builder.py:245-247!
