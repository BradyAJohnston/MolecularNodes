# Node Builder API Examples

## Simple Example

```python
from molecularnodes.nodes import TreeBuilder
from molecularnodes.nodes.sockets import Geometry

tree = TreeBuilder("SimpleOffset")
tree.interface(
    inputs=[Geometry(name="Geometry")],
    outputs=[Geometry(name="Geometry")]
)

with tree:
    tree.input_socket >> SetPosition(offset=(1, 0, 0)) >> tree.output_socket
```

## Complete Example with Multiple Socket Types

```python
from molecularnodes.nodes import TreeBuilder
from molecularnodes.nodes.sockets import Geometry, Boolean, Vector, Float, Integer

tree = TreeBuilder("ComplexTransform")

tree.interface(
    inputs=[
        Geometry(name="Geometry", description="Input geometry to transform"),
        Boolean(name="Apply", default=True, description="Enable transformation"),
        Vector(
            name="Translation",
            default=(0.0, 0.0, 0.0),
            description="World space offset"
        ),
        Float(
            name="Scale",
            default=1.0,
            min_value=0.0,
            max_value=10.0,
            description="Uniform scale multiplier"
        ),
        Integer(
            name="Seed",
            default=0,
            min_value=0,
            description="Random seed"
        ),
    ],
    outputs=[
        Geometry(name="Geometry", description="Transformed geometry"),
        Integer(name="Count", description="Number of vertices processed"),
    ]
)

with tree:
    pos = Position()

    (
        tree.input_socket
        >> SetPosition(position=pos)
        >> TransformGeometry(translation=(0, 0, 1))
        >> tree.output_socket
    )
```

## Real-World Example: Molecular Nodes Style

```python
from molecularnodes.nodes import TreeBuilder
from molecularnodes.nodes.sockets import Geometry, Boolean, Float, Integer, String

tree = TreeBuilder("Color by Attribute")

tree.interface(
    inputs=[
        Geometry(name="Atoms"),
        String(name="Attribute", default="b_factor"),
        Float(name="Min", default=0.0),
        Float(name="Max", default=1.0),
        Boolean(name="Clamp", default=True),
    ],
    outputs=[
        Geometry(name="Atoms"),
    ]
)

with tree:
    # Build the color mapping node tree
    attr = StoreNamedAttribute(name="color_value")
    color_ramp = ColorRamp()

    (
        tree.input_socket
        >> attr
        >> color_ramp
        >> tree.output_socket
    )
```

## Benefits of This API

### 1. Type Safety
```python
# ✅ Type checker catches this
Boolean(name="Selection", default="true")  # Error: Expected bool, got str

# ✅ IDE autocompletes socket types
from molecularnodes.nodes.sockets import Geo...  # IDE suggests Geometry
```

### 2. Discoverable Properties
```python
# IDE shows you what properties are available for each socket type
Vector(
    name="Offset",
    default=(0, 0, 0),      # ← IDE knows this is tuple[float, float, float]
    min_value=(-10, -10, -10),  # ← IDE suggests this property
    max_value=(10, 10, 10),     # ← and this one
    description="..."           # ← and this one
)
```

### 3. Clean, Readable Interface Definitions
```python
# Old way (verbose, error-prone)
tree.tree.interface.new_socket('Geometry', in_out='INPUT', socket_type='NodeSocketGeometry')
tree.tree.interface.new_socket('Selection', in_out='INPUT', socket_type='NodeSocketBool')

# New way (clean, type-safe)
tree.interface(
    inputs=[
        Geometry(name="Geometry"),
        Boolean(name="Selection"),
    ]
)
```

### 4. Grouped and Organized
```python
# Clear separation between inputs and outputs
tree.interface(
    inputs=[
        # All inputs together
        Geometry(name="Geometry"),
        Boolean(name="Selection"),
        Vector(name="Offset"),
    ],
    outputs=[
        # All outputs together
        Geometry(name="Geometry"),
        Integer(name="Count"),
    ]
)
```

### 5. Easy to Generate
```python
# Programmatically create interfaces
input_sockets = [
    Geometry(name="Geometry"),
]

# Add optional parameters based on config
if config.needs_selection:
    input_sockets.append(Boolean(name="Selection", default=True))

if config.needs_offset:
    input_sockets.append(Vector(name="Offset", default=(0, 0, 0)))

tree.interface(inputs=input_sockets, outputs=[Geometry(name="Geometry")])
```
