# Socket vs Node Naming Convention

## The Problem

Without a clear naming convention, there's ambiguity between:
- Socket definition classes (for tree interfaces)
- Node classes (that create actual nodes)

For example:
- `Vector()` could mean a Vector socket or a Vector input node
- `Boolean()` could mean a Boolean socket or a Boolean node

## The Solution: Socket Prefix

All socket definition classes use the `Socket` prefix to clearly distinguish them.

### Socket Definition Classes (for tree interfaces)

```python
from molecularnodes.nodes.sockets import (
    SocketGeometry,
    SocketBoolean,
    SocketVector,
    SocketFloat,
    SocketInt,
    SocketString,
    SocketColor,
    SocketMaterial,
    # ... etc
)

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
```

### Node Classes (create actual nodes in the tree)

```python
from molecularnodes.nodes import Vector, Position, SetPosition

with tree:
    # Vector() creates a Vector INPUT node
    offset = Vector(value=(1, 0, 0))

    # Position() creates a Position INPUT node
    pos = Position()

    # SetPosition() creates a Set Position geometry node
    tree.input_socket >> SetPosition(position=pos, offset=offset)
```

## Clear Examples

### Example 1: Tree Interface Definition

```python
from molecularnodes.nodes import TreeBuilder
from molecularnodes.nodes.sockets import SocketGeometry, SocketVector, SocketFloat

tree = TreeBuilder("OffsetGeometry")

# Socket classes define the interface
tree.interface(
    inputs=[
        SocketGeometry(name="Geometry"),
        SocketVector(name="Translation", default=(0, 0, 0)),
        SocketFloat(name="Scale", default=1.0, min_value=0.0),
    ],
    outputs=[
        SocketGeometry(name="Geometry"),
    ]
)
```

### Example 2: Using Nodes Inside the Tree

```python
from molecularnodes.nodes import Vector, Position, SetPosition, TransformGeometry

with tree:
    # Node classes create nodes
    pos = Position()                    # Creates Position input node
    offset = Vector(value=(1, 0, 0))    # Creates Vector input node

    # Build the node graph
    (
        tree.input_socket
        >> SetPosition(position=pos, offset=offset)
        >> TransformGeometry(translation=(0, 0, 1))
        >> tree.output_socket
    )
```

### Example 3: Mixed Usage

```python
from molecularnodes.nodes import TreeBuilder
from molecularnodes.nodes.sockets import SocketGeometry, SocketBoolean
from molecularnodes.nodes import Position, Vector, SetPosition

tree = TreeBuilder("SmartOffset")

# Socket classes for interface
tree.interface(
    inputs=[
        SocketGeometry(name="Geometry"),
        SocketBoolean(name="Use Custom Offset", default=False),
    ],
    outputs=[
        SocketGeometry(name="Geometry"),
    ]
)

with tree:
    # Node classes for nodes
    pos = Position()
    offset = Vector(value=(2, 0, 0))

    tree.input_socket >> SetPosition(position=pos, offset=offset) >> tree.output_socket
```

## Benefits

1. **No Ambiguity**: `SocketVector` vs `Vector` - clear what each does
2. **Import Clarity**:
   ```python
   from molecularnodes.nodes.sockets import SocketVector  # For interfaces
   from molecularnodes.nodes import Vector                # For nodes
   ```
3. **IDE Autocomplete**: Typing `Socket...` shows all socket types for interfaces
4. **Code Generation**: Easy to generate both socket classes and node classes without naming conflicts
5. **Self-Documenting**: The prefix makes it obvious when you're defining an interface vs building a tree

## Naming Pattern Summary

| Purpose | Prefix | Example | Import From |
|---------|--------|---------|-------------|
| Interface socket definition | `Socket` | `SocketGeometry`, `SocketVector` | `molecularnodes.nodes.sockets` |
| Node in the tree | None | `Vector`, `Position`, `SetPosition` | `molecularnodes.nodes` |
| Tree builder | None | `TreeBuilder` | `molecularnodes.nodes` |

This convention will be maintained when code generation is implemented for both socket classes and node classes.
