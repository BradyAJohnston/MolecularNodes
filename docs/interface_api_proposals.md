# Node Group Interface API Proposals

## Current API (Verbose)
```python
tree = TreeBuilder("MyTree")
tree.tree.interface.new_socket('Geometry', in_out='INPUT', socket_type='NodeSocketGeometry')
tree.tree.interface.new_socket('Geometry', in_out='OUTPUT', socket_type='NodeSocketGeometry')
tree.tree.interface.new_socket('Selection', in_out='INPUT', socket_type='NodeSocketBool')
tree.tree.interface.new_socket('Offset', in_out='INPUT', socket_type='NodeSocketVector')
```

**Problems:**
- Exposes internal Blender API (`tree.tree.interface`)
- Verbose and repetitive
- Socket type strings are not discoverable
- Hard to read at a glance

---

## Proposal 1: Simple Method Calls

Add `add_input()` and `add_output()` methods to TreeBuilder with type inference.

```python
tree = TreeBuilder("MyTree")
tree.add_input("Geometry")  # Type inferred from name
tree.add_input("Selection", "BOOLEAN")
tree.add_input("Offset", "VECTOR")
tree.add_output("Geometry")
```

**Pros:**
- Clean and readable
- Type names are more intuitive than socket type strings
- Easy to understand

**Cons:**
- Still somewhat imperative
- Inputs/outputs scattered across multiple lines
- Type inference by name could be fragile

---

## Proposal 2: Declarative Dictionary/Kwargs

Define inputs and outputs upfront as structured data.

```python
tree = TreeBuilder(
    "MyTree",
    inputs={
        "Geometry": "GEOMETRY",
        "Selection": "BOOLEAN",
        "Offset": "VECTOR"
    },
    outputs={
        "Geometry": "GEOMETRY",
        "Result": "FLOAT"
    }
)
```

**Pros:**
- All interface definition in one place
- Easy to see the full signature at a glance
- Familiar pattern (like function signatures)

**Cons:**
- Dictionaries lose ordering (though Python 3.7+ maintains it)
- Can't set default values or other socket properties easily
- Mixes constructor concerns

---

## Proposal 3: Builder Pattern with Chaining

Use method chaining before entering the tree context.

```python
tree = (
    TreeBuilder("MyTree")
    .input("Geometry")
    .input("Selection", "BOOLEAN")
    .input("Offset", "VECTOR")
    .output("Geometry")
    .output("Result", "FLOAT")
)

with tree:
    # Build nodes...
```

**Pros:**
- Fluent interface, reads like a sentence
- Clear separation between interface and implementation
- Easy to chain

**Cons:**
- Mixes input/output order (though logical grouping)
- Can become verbose for many sockets

---

## Proposal 4: Separate Context Manager for Interface

Use a dedicated context manager for defining the interface.

```python
tree = TreeBuilder("MyTree")

with tree.interface() as iface:
    iface.input("Geometry")
    iface.input("Selection", "BOOLEAN")
    iface.input("Offset", "VECTOR")
    iface.output("Geometry")
    iface.output("Result", "FLOAT")

with tree:
    # Build nodes...
```

**Pros:**
- Clear separation of concerns
- Groups interface definition together
- Could auto-arrange interface in __exit__

**Cons:**
- Two context managers required
- More boilerplate

---

## Proposal 5: Typed Socket Classes (Most Pythonic)

Import typed socket classes and pass them directly.

```python
from molecularnodes.nodes import TreeBuilder, Input, Output
from molecularnodes.nodes.sockets import Geometry, Boolean, Vector, Float

tree = TreeBuilder("MyTree")
tree.interface(
    Input("Geometry", Geometry),
    Input("Selection", Boolean),
    Input("Offset", Vector),
    Output("Geometry", Geometry),
    Output("Result", Float)
)

with tree:
    # Build nodes...
```

**Pros:**
- Fully typed - IDE autocomplete for socket types
- Explicit and clear
- Easy to extend with socket properties
- Type checking works

**Cons:**
- More imports needed
- Slightly more verbose
- Requires defining socket type classes

---

## Proposal 6: Hybrid - Auto Geometry, Explicit Others

Make the common case (single Geometry in/out) automatic, require explicit definition for others.

```python
# Simple case - automatic Geometry in/out
tree = TreeBuilder("MyTree")
# Automatically creates Geometry input/output

# Complex case - explicit interface
tree = TreeBuilder("MyTree", auto_geometry=False)
tree.add_input("Geometry", "GEOMETRY")
tree.add_input("Selection", "BOOLEAN")
tree.add_output("Geometry", "GEOMETRY")
tree.add_output("Count", "INT")
```

**Pros:**
- Minimal boilerplate for common cases
- Still flexible for complex cases
- Progressive disclosure of complexity

**Cons:**
- Magic behavior might be surprising
- Inconsistent API between simple and complex cases

---

## Proposal 7: List of Tuples (Compact)

Use simple data structures for interface definition.

```python
tree = TreeBuilder(
    "MyTree",
    inputs=[
        ("Geometry", "GEOMETRY"),
        ("Selection", "BOOLEAN", False),  # with default value
        ("Offset", "VECTOR", (0, 0, 0))
    ],
    outputs=[
        ("Geometry", "GEOMETRY"),
        ("Result", "FLOAT")
    ]
)
```

**Pros:**
- Compact and readable
- Ordered and explicit
- Easy to include default values
- Familiar tuple-based pattern

**Cons:**
- Tuples are less discoverable than named parameters
- Type checking harder without TypedDict

---

## Recommendation: Proposal 5 (Typed Socket Classes)

This approach offers the best balance of:
- **IDE Support**: Full autocomplete for socket types
- **Type Safety**: Proper type checking
- **Readability**: Clear and explicit
- **Extensibility**: Easy to add socket properties (defaults, min/max, etc.)
- **Pythonic**: Follows Python best practices

### Full Example with Proposal 5:

```python
from molecularnodes.nodes import TreeBuilder, Input, Output
from molecularnodes.nodes.sockets import Geometry, Boolean, Vector

tree = TreeBuilder("MyTree")
tree.interface(
    Input("Geometry", Geometry),
    Input("Selection", Boolean, default=True),
    Input("Offset", Vector, default=(0, 0, 0)),
    Output("Geometry", Geometry)
)

with tree:
    (
        tree.input_socket
        >> SetPosition(offset=tree.input("Offset"))
        >> tree.output_socket
    )
```

This could be further simplified with a helper:

```python
# Even cleaner - socket types inferred from capitalized names
tree.interface(
    Input.Geometry(),
    Input.Selection(default=True),
    Input.Offset(default=(0, 0, 0)),
    Output.Geometry()
)
```

---

## Alternative Recommendation: Proposal 3 (Builder Pattern)

If you prefer less ceremony and fewer imports:

```python
tree = (
    TreeBuilder("MyTree")
    .input("Geometry")
    .input("Selection", "BOOLEAN", default=True)
    .input("Offset", "VECTOR", default=(0, 0, 0))
    .output("Geometry")
)

with tree:
    # Build nodes...
```

This is more concise while still being explicit and readable.
