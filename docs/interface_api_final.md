# Node Group Interface API - Final Design

## API Specification

```python
from molecularnodes.nodes import TreeBuilder
from molecularnodes.nodes.sockets import Geometry, Boolean, Vector, Float, Integer

tree = TreeBuilder("MyTree")
tree.interface(
    inputs=[
        Geometry(name="Geometry"),
        Boolean(name="Selection", default=True),
        Vector(name="Offset", default=(0, 0, 0), min_value=(-10, -10, -10), max_value=(10, 10, 10)),
        Float(name="Scale", default=1.0, min_value=0.0, max_value=5.0),
    ],
    outputs=[
        Geometry(name="Geometry"),
        Integer(name="Count"),
    ]
)

with tree:
    # Access inputs with autocomplete
    tree.input.geometry      # Returns NodeBuilder wrapping Group Input
    tree.input.selection     # Accessing the socket by normalized name
    tree.output.geometry     # Returns NodeBuilder wrapping Group Output
```

## Socket Class Structure

Each socket type is a dataclass with type-specific properties:

```python
from dataclasses import dataclass
from typing import Any

@dataclass
class SocketBase:
    """Base class for all socket definitions."""
    name: str
    bl_socket_type: str = ""
    description: str = ""

@dataclass
class Geometry(SocketBase):
    """Geometry socket - no default value."""
    bl_socket_type: str = "NodeSocketGeometry"

@dataclass
class Boolean(SocketBase):
    """Boolean socket."""
    default: bool = False
    bl_socket_type: str = "NodeSocketBool"

@dataclass
class Vector(SocketBase):
    """Vector socket."""
    default: tuple[float, float, float] = (0.0, 0.0, 0.0)
    min_value: tuple[float, float, float] | None = None
    max_value: tuple[float, float, float] | None = None
    bl_socket_type: str = "NodeSocketVector"

@dataclass
class Float(SocketBase):
    """Float socket."""
    default: float = 0.0
    min_value: float | None = None
    max_value: float | None = None
    bl_socket_type: str = "NodeSocketFloat"

@dataclass
class Integer(SocketBase):
    """Integer socket."""
    default: int = 0
    min_value: int | None = None
    max_value: int | None = None
    bl_socket_type: str = "NodeSocketInt"

@dataclass
class String(SocketBase):
    """String socket."""
    default: str = ""
    bl_socket_type: str = "NodeSocketString"

@dataclass
class Color(SocketBase):
    """Color socket (RGBA)."""
    default: tuple[float, float, float, float] = (1.0, 1.0, 1.0, 1.0)
    bl_socket_type: str = "NodeSocketColor"

@dataclass
class Material(SocketBase):
    """Material socket."""
    bl_socket_type: str = "NodeSocketMaterial"

@dataclass
class Image(SocketBase):
    """Image socket."""
    bl_socket_type: str = "NodeSocketImage"

@dataclass
class Object(SocketBase):
    """Object socket."""
    bl_socket_type: str = "NodeSocketObject"

@dataclass
class Collection(SocketBase):
    """Collection socket."""
    bl_socket_type: str = "NodeSocketCollection"

@dataclass
class Rotation(SocketBase):
    """Rotation socket (Euler or Quaternion)."""
    default: tuple[float, float, float, float] = (1.0, 0.0, 0.0, 0.0)  # Quaternion
    bl_socket_type: str = "NodeSocketRotation"
```

## TreeBuilder Implementation

```python
class TreeBuilder:
    def interface(
        self,
        inputs: list[SocketBase] | None = None,
        outputs: list[SocketBase] | None = None
    ) -> None:
        """Define the node group interface with typed socket definitions.

        Args:
            inputs: List of input socket definitions
            outputs: List of output socket definitions
        """
        if inputs:
            for socket_def in inputs:
                self._create_input_socket(socket_def)

        if outputs:
            for socket_def in outputs:
                self._create_output_socket(socket_def)

    def _create_input_socket(self, socket_def: SocketBase) -> None:
        """Create an input socket from a socket definition."""
        socket = self.tree.interface.new_socket(
            name=socket_def.name,
            in_out='INPUT',
            socket_type=socket_def.bl_socket_type
        )
        self._configure_socket(socket, socket_def)

    def _create_output_socket(self, socket_def: SocketBase) -> None:
        """Create an output socket from a socket definition."""
        socket = self.tree.interface.new_socket(
            name=socket_def.name,
            in_out='OUTPUT',
            socket_type=socket_def.bl_socket_type
        )
        self._configure_socket(socket, socket_def)

    def _configure_socket(self, socket, socket_def: SocketBase) -> None:
        """Configure socket properties from definition."""
        # Set default value if it exists
        if hasattr(socket_def, 'default') and hasattr(socket, 'default_value'):
            socket.default_value = socket_def.default

        # Set min/max values if they exist
        if hasattr(socket_def, 'min_value') and socket_def.min_value is not None:
            if hasattr(socket, 'min_value'):
                socket.min_value = socket_def.min_value

        if hasattr(socket_def, 'max_value') and socket_def.max_value is not None:
            if hasattr(socket, 'max_value'):
                socket.max_value = socket_def.max_value

        # Set description
        if socket_def.description:
            socket.description = socket_def.description
```

## Usage Examples

### Simple Example
```python
from molecularnodes.nodes import TreeBuilder
from molecularnodes.nodes.sockets import Geometry

tree = TreeBuilder("SimpleTree")
tree.interface(
    inputs=[Geometry(name="Geometry")],
    outputs=[Geometry(name="Geometry")]
)

with tree:
    tree.input_socket >> tree.output_socket
```

### Complex Example with Multiple Socket Types
```python
from molecularnodes.nodes import TreeBuilder
from molecularnodes.nodes.sockets import (
    Geometry, Boolean, Vector, Float, Integer, Color, String
)

tree = TreeBuilder("ComplexTree")
tree.interface(
    inputs=[
        Geometry(name="Geometry"),
        Boolean(name="Apply Transform", default=True),
        Vector(
            name="Translation",
            default=(0, 0, 0),
            description="World space translation"
        ),
        Float(
            name="Scale",
            default=1.0,
            min_value=0.0,
            max_value=10.0,
            description="Uniform scale factor"
        ),
        Integer(
            name="Subdivisions",
            default=2,
            min_value=0,
            max_value=6
        ),
        Color(
            name="Tint",
            default=(1.0, 1.0, 1.0, 1.0)
        ),
    ],
    outputs=[
        Geometry(name="Geometry"),
        Integer(name="Vertex Count"),
    ]
)
```

### With Descriptions for Better UX
```python
tree.interface(
    inputs=[
        Geometry(
            name="Geometry",
            description="Input geometry to process"
        ),
        Boolean(
            name="Selection",
            default=True,
            description="Filter which elements to affect"
        ),
    ],
    outputs=[
        Geometry(
            name="Geometry",
            description="Processed geometry"
        ),
    ]
)
```

## Benefits

1. **Full IDE Autocomplete**: Import socket types, IDE suggests them
2. **Type-Safe**: Can't pass wrong type for defaults (e.g., `Boolean(default=1.0)` is a type error)
3. **Self-Documenting**: Socket class names clearly indicate type
4. **Extensible**: Easy to add new socket properties
5. **Code Generation Ready**: Socket classes can be generated from Blender's registry
6. **Clean Separation**: Clear distinction between inputs and outputs
7. **List-Based**: Natural Python pattern, easy to iterate/modify programmatically

## Future Enhancements

1. **Socket Groups**: Support for panel/group organization
2. **Attribute Domains**: For geometry attribute sockets
3. **Custom Properties**: For addon-specific metadata
4. **Validation**: Runtime validation of socket configurations
