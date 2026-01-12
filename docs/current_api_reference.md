# Current Node Builder API Reference

## Complete Example

```python
from molecularnodes.nodes import TreeBuilder
from molecularnodes.nodes.sockets import SocketGeometry, SocketBoolean, SocketVector
from molecularnodes.nodes import Position, Vector, SetPosition, TransformGeometry

# Create a node tree
tree = TreeBuilder("MyNodeTree")

# Define the interface with Socket classes
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

# Build the node graph inside the context manager
with tree:
    # Create nodes (no tree parameter needed - uses active tree)
    pos = Position()
    offset = Vector(value=(1, 0, 0))

    # Chain nodes with >> operator
    # Access specific sockets by name
    (
        tree.input.geometry
        >> SetPosition(
            selection=tree.input.selection,
            position=pos,
            offset=offset
        )
        >> TransformGeometry(translation=(0, 0, 1))
        >> tree.output.geometry
    )
```

## API Components

### 1. TreeBuilder

**Purpose:** Manages a Blender geometry node tree

```python
from molecularnodes.nodes import TreeBuilder

tree = TreeBuilder("TreeName")  # Create new tree
# or
tree = TreeBuilder(existing_tree)  # Wrap existing tree
```

**Context Manager:**
```python
with tree:
    # Nodes created here use the active tree automatically
    node = Position()  # No tree parameter needed
```

**Exit behavior:** Automatically arranges nodes when exiting the context

### 2. Socket Definition Classes

**Purpose:** Define typed interface sockets with properties

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
    SocketImage,
    SocketObject,
    SocketCollection,
    SocketRotation,
    SocketMatrix,
    SocketTexture,
)
```

**Common parameters:**
- `name`: Socket display name (required)
- `default`: Default value (type-specific)
- `min_value`: Minimum value (numeric types)
- `max_value`: Maximum value (numeric types)
- `description`: Tooltip description

**Examples:**
```python
SocketGeometry(name="Geometry")
SocketBoolean(name="Selection", default=True)
SocketVector(name="Offset", default=(0, 0, 0), min_value=(-10, -10, -10))
SocketFloat(name="Scale", default=1.0, min_value=0.0, max_value=5.0)
SocketInt(name="Count", default=10, min_value=0, max_value=100)
```

### 3. TreeBuilder.interface()

**Purpose:** Define tree inputs and outputs

```python
tree.interface(
    inputs=[...],   # List of Socket definition classes
    outputs=[...]   # List of Socket definition classes
)
```

**Example:**
```python
tree.interface(
    inputs=[
        SocketGeometry(name="Geometry", description="Input mesh"),
        SocketBoolean(name="Apply", default=True),
    ],
    outputs=[
        SocketGeometry(name="Geometry", description="Output mesh"),
    ]
)
```

### 4. Node Classes

**Purpose:** Create actual nodes in the tree

Currently available:
- `Position()` - Input Position node
- `Vector(value=(x, y, z))` - Input Vector node
- `SetPosition(...)` - Set Position geometry node
- `TransformGeometry(...)` - Transform Geometry node

**Usage:**
```python
with tree:
    # No tree parameter - uses active context
    pos = Position()
    vec = Vector(value=(1, 0, 0))

    # With inputs
    set_pos = SetPosition(
        selection=some_boolean,
        position=pos,
        offset=vec
    )
```

### 5. Node Chaining with >> Operator

**Purpose:** Connect nodes together

```python
node1 >> node2 >> node3
```

The `>>` operator:
- Prefers "Geometry" sockets for geometry nodes
- Falls back to default input/output if no Geometry socket
- Returns the right-hand node for continued chaining

**Example:**
```python
with tree:
    (
        tree.input_socket
        >> SetPosition(position=Position())
        >> TransformGeometry(translation=(0, 0, 1))
        >> tree.output_socket
    )
```

### 6. Named Socket Access

**Access tree sockets by name:**
- `tree.input.socket_name` - Access input socket by normalized name
- `tree.output.socket_name` - Access output socket by normalized name

**Name normalization:**
Socket names are normalized to Python attribute names:
- `"Geometry"` → `.geometry`
- `"Selection"` → `.selection`
- `"My Socket"` → `.my_socket`
- `"Point Count"` → `.point_count`

**Examples:**
```python
# Single socket
tree.interface(inputs=[SocketGeometry(name="Geometry")])
with tree:
    tree.input.geometry >> some_node

# Multiple sockets
tree.interface(inputs=[
    SocketGeometry(name="Geometry"),
    SocketBoolean(name="Selection"),
])
with tree:
    tree.input.geometry >> SetPosition(selection=tree.input.selection)
```


## Naming Conventions

### Socket Classes vs Node Classes

**Socket classes** (for interface definitions):
- Prefix: `Socket`
- Import from: `molecularnodes.nodes.sockets`
- Example: `SocketVector`, `SocketGeometry`
- Purpose: Define interface socket properties

**Node classes** (for nodes in the tree):
- Prefix: None
- Import from: `molecularnodes.nodes`
- Example: `Vector`, `Position`, `SetPosition`
- Purpose: Create actual nodes in the tree

### Why the distinction?

```python
# Clear and unambiguous:
from molecularnodes.nodes.sockets import SocketVector  # Interface socket
from molecularnodes.nodes import Vector                # Vector input node

tree.interface(inputs=[SocketVector(name="Offset")])   # Define interface
with tree:
    offset = Vector(value=(1, 0, 0))                   # Create node
```

## What's Next

### Pending features:
1. Named access to tree inputs/outputs (e.g., `tree.inputs.geometry`)
2. Code generator for creating node classes from Blender registry
3. Socket wrapper classes with type hints
4. More generated node classes (currently only 4 manual examples)

### Current limitations:
- Limited node classes (only Position, Vector, SetPosition, TransformGeometry)
- Generic input/output socket access (no per-socket access yet)
- Manual node class definitions (no code generation yet)
