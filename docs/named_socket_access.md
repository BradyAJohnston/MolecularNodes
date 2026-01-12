# Named Socket Access

## The Improvement

We've updated socket access to be more explicit and discoverable.

### Before (Generic Access)
```python
tree.interface(
    inputs=[SocketGeometry(name="Geometry")],
    outputs=[SocketGeometry(name="Geometry")]
)

with tree:
    tree.input_socket >> SetPosition() >> tree.output_socket
    # ❌ Which socket? Not clear when you have multiple
```

### After (Named Access)
```python
tree.interface(
    inputs=[SocketGeometry(name="Geometry")],
    outputs=[SocketGeometry(name="Geometry")]
)

with tree:
    tree.input.geometry >> SetPosition() >> tree.output.geometry
    # ✅ Explicit about which socket you're accessing
```

## How It Works

### Name Normalization

Socket names are automatically normalized for Python attribute access:

| Socket Name | Python Attribute |
|-------------|------------------|
| `"Geometry"` | `.geometry` |
| `"Selection"` | `.selection` |
| `"My Socket"` | `.my_socket` |
| `"Point Count"` | `.point_count` |

### SocketAccessor Class

`TreeBuilder` has two `SocketAccessor` instances:
- `tree.input` - Access input sockets
- `tree.output` - Access output sockets

When you access an attribute like `tree.input.geometry`:
1. The name is denormalized: `"geometry"` → `"Geometry"`
2. A `SocketNodeBuilder` is created wrapping that specific socket
3. The `SocketNodeBuilder` works with the `>>` operator

## Examples

### Simple Case: Single Geometry Socket

```python
from molecularnodes.nodes import TreeBuilder
from molecularnodes.nodes.sockets import SocketGeometry

tree = TreeBuilder("SimpleTree")
tree.interface(
    inputs=[SocketGeometry(name="Geometry")],
    outputs=[SocketGeometry(name="Geometry")]
)

with tree:
    tree.input.geometry >> SetPosition() >> tree.output.geometry
```

### Multiple Inputs

```python
from molecularnodes.nodes.sockets import SocketGeometry, SocketBoolean, SocketVector

tree = TreeBuilder("MultiInputTree")
tree.interface(
    inputs=[
        SocketGeometry(name="Geometry"),
        SocketBoolean(name="Selection", default=True),
        SocketVector(name="Offset", default=(0, 0, 0)),
    ],
    outputs=[
        SocketGeometry(name="Geometry"),
    ]
)

with tree:
    # Access each input by name
    (
        tree.input.geometry
        >> SetPosition(
            selection=tree.input.selection,
            offset=tree.input.offset
        )
        >> tree.output.geometry
    )
```

### Multiple Outputs

```python
from molecularnodes.nodes.sockets import SocketGeometry, SocketInt

tree = TreeBuilder("MultiOutputTree")
tree.interface(
    inputs=[SocketGeometry(name="Geometry")],
    outputs=[
        SocketGeometry(name="Geometry"),
        SocketInt(name="Vertex Count"),
        SocketInt(name="Face Count"),
    ]
)

with tree:
    # Connect to different output sockets
    tree.input.geometry >> some_node >> tree.output.geometry

    # Could connect other nodes to:
    # tree.output.vertex_count
    # tree.output.face_count
```

### Custom Socket Names

```python
tree.interface(
    inputs=[
        SocketGeometry(name="Base Mesh"),
        SocketGeometry(name="Deform Cage"),
    ],
    outputs=[
        SocketGeometry(name="Result"),
    ]
)

with tree:
    # Normalized names with underscores
    base = tree.input.base_mesh
    cage = tree.input.deform_cage

    base >> some_node >> tree.output.result
```

## Benefits

1. **Explicit**: Clear which socket you're connecting to
2. **Discoverable**: Type `tree.input.` and see available sockets (via __getattr__)
3. **Scalable**: Works with any number of inputs/outputs
4. **Type-safe ready**: Foundation for adding type hints later
5. **Consistent**: Matches the declarative interface definition style

## Implementation Details

### SocketAccessor
```python
class SocketAccessor:
    def __init__(self, tree: TreeBuilder, direction: str):
        self._tree = tree
        self._direction = direction  # 'INPUT' or 'OUTPUT'

    def __getattr__(self, name: str) -> NodeBuilder:
        socket_name = denormalize_name(name)  # geometry -> Geometry
        node = self._tree.input() if self._direction == "INPUT" else self._tree.output()
        return SocketNodeBuilder(node, socket_name, self._direction)
```

### SocketNodeBuilder
```python
class SocketNodeBuilder(NodeBuilder):
    """Wraps a specific socket on an input/output node."""

    @property
    def default_output(self) -> NodeSocket:
        if self._direction == "INPUT":
            return self.node.outputs[self._socket_name]
        else:
            return self.node.outputs[0]

    @property
    def default_input(self) -> NodeSocket:
        if self._direction == "OUTPUT":
            return self.node.inputs[self._socket_name]
        else:
            return self.node.inputs[0]
```

This ensures the `>>` operator connects to the right socket.
