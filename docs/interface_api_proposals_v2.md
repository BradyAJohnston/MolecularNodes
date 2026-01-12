# Node Group Interface API Proposals v2
## Design Principles: Type Safety & IDE Autocomplete

All proposals avoid string-based type specification and provide full IDE autocomplete.

---

## Proposal 1: Typed Input/Output Classes

Define socket types as classes that can be imported and used directly.

```python
from molecularnodes.nodes import TreeBuilder
from molecularnodes.nodes.sockets import Geometry, Boolean, Vector, Float, Integer

tree = TreeBuilder("MyTree")

# Option 1a: Class-based with kwargs
tree.add_input(Geometry(name="Geometry"))
tree.add_input(Boolean(name="Selection", default=True))
tree.add_input(Vector(name="Offset", default=(0, 0, 0)))
tree.add_output(Geometry(name="Geometry"))

# Option 1b: Using Input/Output wrapper classes
from molecularnodes.nodes.interface import Input, Output

tree.interface(
    Input.geometry("Geometry"),
    Input.boolean("Selection", default=True),
    Input.vector("Offset", default=(0, 0, 0)),
    Output.geometry("Geometry")
)
```

**Implementation:**
```python
# molecularnodes/nodes/interface.py
class SocketType:
    def __init__(self, name: str, default=None):
        self.name = name
        self.default = default

class Geometry(SocketType):
    bl_socket_type = "NodeSocketGeometry"

class Boolean(SocketType):
    bl_socket_type = "NodeSocketBool"
    def __init__(self, name: str, default: bool = False):
        super().__init__(name, default)

class Input:
    @staticmethod
    def geometry(name: str) -> Geometry:
        return Geometry(name)

    @staticmethod
    def boolean(name: str, default: bool = False) -> Boolean:
        return Boolean(name, default)
```

**Pros:**
- Full type hints and IDE autocomplete
- No string-based type specification
- Default values are type-checked
- Extensible for socket properties

**Cons:**
- Requires defining socket type classes
- More verbose than string-based approach

---

## Proposal 2: Generic Input/Output with Type Parameters

Use Python's typing system with generic classes.

```python
from molecularnodes.nodes import TreeBuilder, Input, Output
from molecularnodes.nodes.types import Geometry, Boolean, Vector

tree = TreeBuilder("MyTree")
tree.interface(
    Input[Geometry]("Geometry"),
    Input[Boolean]("Selection", default=True),
    Input[Vector]("Offset", default=(0, 0, 0)),
    Output[Geometry]("Geometry")
)

with tree:
    # Access inputs with autocomplete
    sel = tree.input.selection  # Type: Boolean
    offset = tree.input.offset  # Type: Vector
```

**Implementation:**
```python
from typing import Generic, TypeVar

T = TypeVar('T', bound=SocketType)

class Input(Generic[T]):
    def __init__(self, name: str, default: T | None = None):
        self.name = name
        self.default = default
        self.socket_type: type[T] = self.__orig_class__.__args__[0]
```

**Pros:**
- Leverages Python's type system
- Full IDE autocomplete
- Type checking works properly
- Clean, Pythonic syntax

**Cons:**
- Generic syntax might be unfamiliar
- Requires Python 3.9+ for clean syntax
- Runtime type extraction is tricky

---

## Proposal 3: Builder Pattern with Typed Methods

Use method chaining with strongly typed methods.

```python
from molecularnodes.nodes import TreeBuilder

tree = (
    TreeBuilder("MyTree")
    .with_geometry_input(name="Geometry")
    .with_boolean_input(name="Selection", default=True)
    .with_vector_input(name="Offset", default=(0, 0, 0))
    .with_geometry_output(name="Geometry")
)

with tree:
    # Access with autocomplete
    tree.input.geometry  # Returns GeometrySocket
    tree.input.selection  # Returns BooleanSocket
```

**Implementation:**
```python
class TreeBuilder:
    def with_geometry_input(self, name: str) -> "TreeBuilder":
        self.tree.interface.new_socket(name, in_out='INPUT', socket_type='NodeSocketGeometry')
        return self

    def with_boolean_input(self, name: str, default: bool = False) -> "TreeBuilder":
        socket = self.tree.interface.new_socket(name, in_out='INPUT', socket_type='NodeSocketBool')
        socket.default_value = default
        return self
```

**Pros:**
- Method names provide type information
- Full autocomplete on method names
- Type hints on parameters
- Familiar builder pattern

**Cons:**
- Requires one method per socket type
- Method names can be verbose (with_X_input)
- Doesn't scale to hundreds of socket types

---

## Proposal 4: Enum-Based Socket Types

Use Python enums for socket types with factory methods.

```python
from molecularnodes.nodes import TreeBuilder, Input, Output, SocketType

tree = TreeBuilder("MyTree")
tree.interface(
    Input(SocketType.GEOMETRY, "Geometry"),
    Input(SocketType.BOOLEAN, "Selection", default=True),
    Input(SocketType.VECTOR, "Offset", default=(0, 0, 0)),
    Output(SocketType.GEOMETRY, "Geometry")
)
```

**Implementation:**
```python
from enum import Enum

class SocketType(Enum):
    GEOMETRY = "NodeSocketGeometry"
    BOOLEAN = "NodeSocketBool"
    VECTOR = "NodeSocketVector"
    FLOAT = "NodeSocketFloat"
    INTEGER = "NodeSocketInt"
    # ... etc

class Input:
    def __init__(self, socket_type: SocketType, name: str, default=None):
        self.socket_type = socket_type
        self.name = name
        self.default = default
```

**Pros:**
- Enums provide autocomplete
- Type-safe socket type specification
- Standard Python pattern
- Easy to add new types

**Cons:**
- Still somewhat string-like (enum values)
- Doesn't provide type-specific parameter validation
- Less granular type hints

---

## Proposal 5: Protocol-Based with Dataclasses

Use dataclasses for socket definitions with protocols for type checking.

```python
from molecularnodes.nodes import TreeBuilder
from molecularnodes.nodes.sockets import GeometryInput, BooleanInput, VectorInput, GeometryOutput

tree = TreeBuilder("MyTree")
tree.interface(
    GeometryInput(name="Geometry"),
    BooleanInput(name="Selection", default=True),
    VectorInput(name="Offset", default=(0, 0, 0)),
    GeometryOutput(name="Geometry")
)

with tree:
    # Strongly typed access
    tree.inputs.geometry  # Type: GeometrySocket
    tree.inputs.selection  # Type: BoolSocket
```

**Implementation:**
```python
from dataclasses import dataclass
from typing import Protocol

class SocketDefinition(Protocol):
    name: str
    bl_socket_type: str

@dataclass
class GeometryInput:
    name: str
    bl_socket_type: str = "NodeSocketGeometry"

@dataclass
class BooleanInput:
    name: str
    default: bool = False
    bl_socket_type: str = "NodeSocketBool"

@dataclass
class VectorInput:
    name: str
    default: tuple[float, float, float] = (0.0, 0.0, 0.0)
    bl_socket_type: str = "NodeSocketVector"
```

**Pros:**
- Dataclasses provide great IDE support
- Each socket type has its own class
- Type-specific default values
- Clean, declarative syntax
- Can validate defaults at init time

**Cons:**
- Many socket classes to define (but could be generated!)
- Separate classes for Input vs Output

---

## Proposal 6: Functional API with Type Classes

Use simple functions that return typed objects.

```python
from molecularnodes.nodes import TreeBuilder
from molecularnodes.nodes import sockets

tree = TreeBuilder("MyTree")
tree.interface(
    sockets.geometry_input("Geometry"),
    sockets.boolean_input("Selection", default=True),
    sockets.vector_input("Offset", default=(0, 0, 0)),
    sockets.geometry_output("Geometry")
)
```

**Implementation:**
```python
# molecularnodes/nodes/sockets.py
@dataclass
class GeometryInput:
    name: str
    bl_socket_type: str = "NodeSocketGeometry"
    direction: str = "INPUT"

def geometry_input(name: str) -> GeometryInput:
    """Create a geometry input socket."""
    return GeometryInput(name=name)

def boolean_input(name: str, default: bool = False) -> BooleanInput:
    """Create a boolean input socket."""
    return BooleanInput(name=name, default=default)
```

**Pros:**
- Functions show up in autocomplete
- Docstrings visible in IDE
- Type hints on return values
- Clean, functional style

**Cons:**
- More functions to define
- Namespace could get crowded

---

## Proposal 7: Class Properties with Decorator

Most innovative - use a decorator to define interface declaratively.

```python
from molecularnodes.nodes import node_tree
from molecularnodes.nodes.sockets import Geometry, Boolean, Vector

@node_tree("MyTree")
class MyTree:
    # Inputs
    geometry: Geometry
    selection: Boolean = True
    offset: Vector = (0, 0, 0)

    # Outputs
    out_geometry: Geometry

    def build(self, tree: TreeBuilder):
        with tree:
            (
                tree.input.geometry
                >> SetPosition(offset=tree.input.offset)
                >> tree.output.out_geometry
            )
```

**Implementation:**
```python
def node_tree(name: str):
    def decorator(cls):
        # Introspect annotations to create interface
        tree = TreeBuilder(name)
        for field_name, field_type in cls.__annotations__.items():
            if field_name.startswith('out_'):
                # Output
                default = getattr(cls, field_name, None)
                tree._add_output(field_type, field_name, default)
            else:
                # Input
                default = getattr(cls, field_name, None)
                tree._add_input(field_type, field_name, default)

        instance = cls()
        instance.build(tree)
        return tree
    return decorator
```

**Pros:**
- Type annotations provide full IDE support
- Declarative and clean
- Defaults work naturally
- Very Pythonic (like dataclasses)

**Cons:**
- Magic/unusual pattern for Blender
- Harder to understand for beginners
- Mixing class definition with node tree building

---

## Recommendation: Proposal 5 (Dataclass-Based)

This provides the best balance:

```python
from molecularnodes.nodes import TreeBuilder
from molecularnodes.nodes.sockets import (
    GeometryInput,
    BooleanInput,
    VectorInput,
    GeometryOutput
)

tree = TreeBuilder("MyTree")
tree.interface(
    GeometryInput(name="Geometry"),
    BooleanInput(name="Selection", default=True),
    VectorInput(name="Offset", default=(0, 0, 0)),
    GeometryOutput(name="Geometry")
)

with tree:
    # Access with full type hints
    (
        tree.input.geometry
        >> SetPosition(
            selection=tree.input.selection,
            offset=tree.input.offset
        )
        >> tree.output.geometry
    )
```

### Why This Approach?

1. **Full IDE autocomplete**: Import socket types, IDE suggests them
2. **Type-specific validation**: `BooleanInput(default=True)` - can't pass wrong type
3. **Generated code friendly**: Socket classes can be code-generated just like nodes
4. **Clean syntax**: Reads naturally
5. **Extensible**: Easy to add socket properties (min/max, descriptions, etc.)

### Implementation Strategy

Generate socket definition classes from Blender's socket types:

```python
# molecularnodes/nodes/sockets.py (generated)
from dataclasses import dataclass

@dataclass
class GeometryInput:
    """Geometry socket input."""
    name: str
    bl_socket_type: str = "NodeSocketGeometry"

@dataclass
class BooleanInput:
    """Boolean socket input."""
    name: str
    default: bool = False
    bl_socket_type: str = "NodeSocketBool"

@dataclass
class VectorInput:
    """Vector socket input."""
    name: str
    default: tuple[float, float, float] = (0.0, 0.0, 0.0)
    bl_socket_type: str = "NodeSocketVector"

# ... etc for all socket types
```

This aligns perfectly with your code generation strategy!
