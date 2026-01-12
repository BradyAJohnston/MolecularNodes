"""
Code generator for Blender geometry node classes.

This script introspects Blender's node registry and generates Python classes
for all geometry nodes with proper type hints and autocomplete support.

Run this script from within Blender to generate node classes:
    blender --background --python generator.py

Or import and run from Blender's Python console:
    import molecularnodes.nodes.generator as gen
    gen.generate_all()
"""

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import bpy


@dataclass
class SocketInfo:
    """Information about a node socket."""

    name: str
    identifier: str  # Internal identifier
    label: str  # Socket label (empty string if no label)
    bl_socket_type: str  # e.g., "NodeSocketGeometry", "NodeSocketFloat"
    socket_type: str  # e.g., "GEOMETRY", "FLOAT", "VECTOR"
    is_output: bool
    is_multi_input: bool = False
    default_value: Any = None
    min_value: Any = None
    max_value: Any = None


@dataclass
class PropertyInfo:
    """Information about a node property."""

    identifier: str
    name: str
    prop_type: str  # 'ENUM', 'BOOLEAN', 'INT', 'FLOAT', etc.
    enum_items: list[tuple[str, str]] | None = None  # [(identifier, name), ...]
    default: Any = None


@dataclass
class NodeInfo:
    """Complete information about a node type."""

    bl_idname: str  # e.g., "GeometryNodeSetPosition"
    name: str  # Human-readable name
    category: str  # e.g., "Geometry", "Mesh", "Curve"
    description: str
    inputs: list[SocketInfo]
    outputs: list[SocketInfo]
    properties: list[PropertyInfo]


def normalize_name(name: str) -> str:
    """Convert 'Geometry' or 'My Socket' to 'geometry' or 'my_socket'.

    Handles numeric names by prefixing with 'input_' to make valid Python identifiers.
    """
    normalized = name.lower().replace(" ", "_").replace("-", "_")

    # If the name starts with a digit or is purely numeric, prefix it
    if normalized and (normalized[0].isdigit() or normalized.isdigit()):
        normalized = f"input_{normalized}"

    return normalized


def get_socket_param_name(socket: SocketInfo, use_identifier: bool = False) -> str:
    """Get the best parameter name for a socket, preferring label over name."""
    # Use label if available and non-empty, otherwise fallback to name
    # if sockets all use the same label name, we need to drop back to using the iden
    if use_identifier:
        return normalize_name(socket.identifier)
    else:
        display_name = socket.label if socket.label else socket.name
        return normalize_name(display_name)


def python_class_name(name: str) -> str:
    """Convert 'set position' to 'SetPosition'.

    Handles edge cases like '3D Cursor' -> 'Cursor3D' to ensure valid Python identifiers.
    """
    # Replace common separators with spaces
    class_name = name.replace("_", " ").replace("-", " ")

    # Remove any characters that aren't alphanumeric or spaces
    # This handles cases like "Field Min&Max" -> "Field Min Max"
    class_name = "".join(c if c.isalnum() or c.isspace() else "" for c in class_name)

    # Title case and remove spaces
    class_name = class_name.title().replace(" ", "")

    # Python class names can't start with a digit, so move leading numbers to the end
    if class_name and class_name[0].isdigit():
        # Find where digits end
        i = 0
        while i < len(class_name) and class_name[i].isdigit():
            i += 1
        # Move digits to end: "3DCursor" -> "Cursor3D"
        class_name = class_name[i:] + class_name[:i]

    return class_name


def get_socket_type_hint(socket_info: SocketInfo) -> str:
    """Get the Python type hint for a socket."""
    # Map Blender socket types to our type hints
    type_map = {
        "NodeSocketGeometry": "LINKABLE",
        "NodeSocketBool": "TYPE_INPUT_BOOLEAN",
        "NodeSocketVector": "TYPE_INPUT_VECTOR",
        "NodeSocketRotation": "TYPE_INPUT_ROTATION",
        "NodeSocketFloat": "float | LINKABLE | None",
        "NodeSocketInt": "int | LINKABLE | None",
        "NodeSocketString": "str | LINKABLE | None",
        "NodeSocketColor": "tuple[float, float, float, float] | LINKABLE | None",
        "NodeSocketMaterial": "LINKABLE | None",
        "NodeSocketImage": "LINKABLE | None",
        "NodeSocketObject": "LINKABLE | None",
        "NodeSocketCollection": "LINKABLE | None",
    }

    return type_map.get(socket_info.bl_socket_type, "LINKABLE | None")


def introspect_node(node_type: type) -> NodeInfo | None:
    """Introspect a Blender node type and extract all information."""
    try:
        # Create temporary node group to instantiate the node
        temp_tree = bpy.data.node_groups.new("temp", "GeometryNodeTree")
        node = temp_tree.nodes.new(node_type.__name__)

        # Extract basic info
        bl_idname = node_type.__name__
        name = node_type.bl_rna.name
        description = node_type.bl_rna.description or f"{name} node"

        # Try to determine category from the node's menu registration
        # This is a simplified approach - geometry-script has more sophisticated categorization
        category = "Geometry"  # Default
        if "Mesh" in bl_idname:
            category = "Mesh"
        elif "Curve" in bl_idname:
            category = "Curve"
        elif "Attribute" in bl_idname:
            category = "Attribute"
        elif "Input" in bl_idname:
            category = "Input"
        elif "Utility" in bl_idname or "Math" in bl_idname:
            category = "Utilities"
        elif bl_idname.startswith("ShaderNode"):
            # Shader nodes used in geometry contexts
            if bl_idname in [
                "ShaderNodeMath",
                "ShaderNodeVectorMath",
                "ShaderNodeClamp",
                "ShaderNodeMapRange",
                "ShaderNodeMix",
                "ShaderNodeInvert",
            ]:
                category = "Utilities"
            elif bl_idname in [
                "ShaderNodeCombineXYZ",
                "ShaderNodeSeparateXYZ",
                "ShaderNodeCombineColor",
                "ShaderNodeSeparateColor",
            ]:
                category = "Utilities"
            elif bl_idname in ["ShaderNodeRGB", "ShaderNodeValue"]:
                category = "Input"
            else:
                category = "Utilities"  # Default for other shader nodes
        elif bl_idname.startswith("FunctionNode"):
            # Function nodes
            if "Input" in bl_idname:
                category = "Input"
            else:
                category = "Utilities"

        # Extract input sockets
        inputs = []
        for socket in node.inputs:
            if not socket.enabled:
                continue

            socket_info = SocketInfo(
                name=socket.name,
                identifier=socket.identifier,
                label=getattr(socket, "label", ""),  # Capture socket label
                bl_socket_type=type(socket).__name__,
                socket_type=socket.type,
                is_output=False,
                is_multi_input=getattr(socket, "is_multi_input", False),
            )

            # Try to get default value
            if hasattr(socket, "default_value"):
                try:
                    socket_info.default_value = socket.default_value
                except (AttributeError, TypeError):
                    pass

            # Try to get min/max
            if hasattr(socket, "min_value"):
                socket_info.min_value = socket.min_value
            if hasattr(socket, "max_value"):
                socket_info.max_value = socket.max_value

            inputs.append(socket_info)

        # Extract output sockets
        outputs = []
        for socket in node.outputs:
            if not socket.enabled:
                continue

            socket_info = SocketInfo(
                name=socket.name,
                identifier=socket.identifier,
                label=getattr(socket, "label", ""),  # Capture socket label
                bl_socket_type=type(socket).__name__,
                socket_type=socket.type,
                is_output=True,
            )
            outputs.append(socket_info)

        # Extract properties (enum menus, boolean flags, etc.)
        properties = []
        parent_props = set()
        for base in node_type.__bases__:
            if hasattr(base, "bl_rna"):
                for prop in base.bl_rna.properties:
                    parent_props.add(prop.identifier)

        for prop in node_type.bl_rna.properties:
            if prop.identifier in parent_props:
                continue

            if prop.type == "ENUM":
                enum_items = [(item.identifier, item.name) for item in prop.enum_items]
                properties.append(
                    PropertyInfo(
                        identifier=prop.identifier,
                        name=prop.name,
                        prop_type="ENUM",
                        enum_items=enum_items,
                    )
                )
            elif prop.type in ["BOOLEAN", "INT", "FLOAT", "STRING"]:
                properties.append(
                    PropertyInfo(
                        identifier=prop.identifier,
                        name=prop.name,
                        prop_type=prop.type,
                    )
                )

        # Clean up
        bpy.data.node_groups.remove(temp_tree)

        return NodeInfo(
            bl_idname=bl_idname,
            name=name,
            category=category,
            description=description,
            inputs=inputs,
            outputs=outputs,
            properties=properties,
        )

    except Exception as e:
        print(f"Error introspecting {node_type.__name__}: {e}")
        return None


def generate_enum_class_methods(node_info: NodeInfo) -> str:
    """Generate @classmethod convenience methods for enum operations."""
    methods = []

    # Find the main operation enum (usually contains "operation" in name)
    operation_enum = None
    for prop in node_info.properties:
        if prop.prop_type == "ENUM" and (
            "operation" in prop.identifier.lower() or prop.identifier == "type"
        ):
            operation_enum = prop
            break

    if not operation_enum:
        return ""

    class_name = python_class_name(node_info.name)

    # Generate method for each enum value
    for enum_id, enum_name in operation_enum.enum_items:
        method_name = enum_id.lower()

        # Add underscore suffix to avoid Python keyword conflicts
        method_name = f"{method_name}_"

        # Skip invalid method names
        if not method_name.replace("_", "").isalnum():
            continue

        # Generate method signature based on node inputs (excluding operation socket)
        input_params = ["cls"]
        call_params = []

        all_labels = [socket.label for socket in node_info.inputs]
        sockets_use_same_name = all(label == all_labels[0] for label in all_labels)
        for socket in node_info.inputs:
            # Use label-based parameter naming
            param_name = get_socket_param_name(socket, sockets_use_same_name)
            if (
                param_name
                and param_name != ""
                and param_name != normalize_name(operation_enum.identifier)
            ):
                type_hint = get_socket_type_hint(socket)
                input_params.append(f"{param_name}: {type_hint} = None")
                call_params.append(f"{param_name}={param_name}")

        params_str = ", ".join(input_params)
        call_params_str = ", ".join(call_params)

        # Add operation parameter to call
        operation_param = f'{normalize_name(operation_enum.identifier)}="{enum_id}"'
        if call_params_str:
            call_params_str = f"{operation_param}, {call_params_str}"
        else:
            call_params_str = operation_param

        method = f'''
    @classmethod
    def {method_name}(
        {params_str}
    ) -> "{class_name}":
        """Create {node_info.name} with operation '{enum_name}'."""
        return cls({call_params_str})'''

        methods.append(method)

    return "".join(methods)


def generate_node_class(node_info: NodeInfo) -> tuple[str, bool]:
    """Generate Python class code for a node.

    Returns:
        tuple[str, bool]: (generated code, has_dynamic_sockets)
    """
    class_name = python_class_name(node_info.name)
    has_dynamic_sockets = False

    # Build __init__ parameters
    init_params = ["self"]
    establish_links_params = []

    # Add input sockets as parameters
    all_labels = [socket.label for socket in node_info.inputs]
    sockets_use_same_name = all(label == all_labels[0] for label in all_labels)
    for socket in node_info.inputs:
        param_name = get_socket_param_name(socket, sockets_use_same_name)

        # Skip unnamed sockets (dynamic sockets that users can drag into)
        # TODO: Support dynamic multi-input sockets properly
        if not param_name or param_name.strip() == "":
            has_dynamic_sockets = True
            continue

        type_hint = get_socket_type_hint(socket)
        init_params.append(f"{param_name}: {type_hint} = None")
        establish_links_params.append((param_name, socket))

    # Add properties as parameters
    for prop in node_info.properties:
        param_name = normalize_name(prop.identifier)
        if prop.prop_type == "ENUM":
            init_params.append(f"{param_name}: str | None = None")
        elif prop.prop_type == "BOOLEAN":
            init_params.append(f"{param_name}: bool | None = None")
        elif prop.prop_type == "INT":
            init_params.append(f"{param_name}: int | None = None")
        elif prop.prop_type == "FLOAT":
            init_params.append(f"{param_name}: float | None = None")
        else:
            init_params.append(f"{param_name}: Any | None = None")

    # Format init signature
    if len(init_params) > 2:  # If more than just self
        init_signature = "(\n        " + ",\n        ".join(init_params) + "\n    )"
    else:
        init_signature = "(" + ", ".join(init_params) + ")"

    # Build establish_links call - map parameter names to socket identifiers
    establish_call = ""
    if establish_links_params:
        # Create mapping of parameter names to socket identifiers for _establish_links
        link_mappings = []
        for param_name, socket in establish_links_params:
            # Use socket identifier as key (which maps to the actual blender socket)
            # parameter name as value
            link_mappings.append(f'"{socket.identifier}": {param_name}')

        if link_mappings:
            establish_call = f"""        self._establish_links(**{{
            {", ".join(link_mappings)}
        }})"""

    # Build property setting calls
    property_calls = []
    for prop in node_info.properties:
        param_name = normalize_name(prop.identifier)
        property_calls.append(f"""        if {param_name} is not None:
            self.node.{prop.identifier} = {param_name}""")

    property_setting = "\n".join(property_calls) if property_calls else ""

    # Generate output properties
    output_properties = []
    used_output_names = set()
    for socket in node_info.outputs:
        # Use socket identifier if name would conflict, otherwise use name
        prop_name = normalize_name(socket.name)
        if prop_name in used_output_names:
            prop_name = (
                normalize_name(socket.identifier)
                if socket.identifier
                else f"{prop_name}_{socket.identifier}"
            )

        # Skip unnamed sockets (dynamic output sockets)
        # TODO: Support dynamic multi-output sockets properly
        if not prop_name or prop_name.strip() == "":
            continue

        used_output_names.add(prop_name)
        output_properties.append(f'''
    @property
    def {prop_name}(self) -> NodeSocket:
        """Output socket: {socket.name}"""
        return self.node.outputs["{socket.name}"]''')

    # Generate enum convenience methods
    enum_methods = generate_enum_class_methods(node_info)

    # Build class
    class_code = f'''
class {class_name}(NodeBuilder):
    """{node_info.description}"""
    name = "{node_info.bl_idname}"

    def __init__{init_signature}:
        super().__init__()
{establish_call}
{property_setting}
{enum_methods}
{"".join(output_properties)}
'''

    return class_code.strip(), has_dynamic_sockets


def get_all_geometry_nodes() -> list[type]:
    """Get all registered geometry node types from Blender, including useful shader and function nodes."""

    # Useful shader nodes that work well in geometry node trees
    USEFUL_SHADER_NODES = [
        "ShaderNodeMath",
        "ShaderNodeVectorMath",
        "ShaderNodeMix",
        "ShaderNodeClamp",
        "ShaderNodeMapRange",
        "ShaderNodeCombineXYZ",
        "ShaderNodeSeparateXYZ",
        "ShaderNodeRGB",
        "ShaderNodeValue",
        "ShaderNodeCombineColor",
        "ShaderNodeSeparateColor",
        "ShaderNodeInvert",
        "ShaderNodeGamma",
        "ShaderNodeBrightContrast",
        "ShaderNodeHueSaturation",
        "ShaderNodeRGBToBW",
        "ShaderNodeFloatCurve",
        "ShaderNodeRGBCurve",
    ]

    # Useful function nodes for utilities
    USEFUL_FUNCTION_NODES = [
        "FunctionNodeCompare",
        "FunctionNodeBooleanMath",
        "FunctionNodeFloatToInt",
        "FunctionNodeRandomValue",
        "FunctionNodeCombineColor",
        "FunctionNodeSeparateColor",
        "FunctionNodeHashValue",
        "FunctionNodeInputBool",
        "FunctionNodeInputColor",
        "FunctionNodeInputInt",
        "FunctionNodeInputString",
        "FunctionNodeInputVector",
    ]

    all_nodes = []

    for attr_name in dir(bpy.types):
        node_type = getattr(bpy.types, attr_name)

        # Check if it's a registered node type
        if not isinstance(node_type, type):
            continue

        if not issubclass(node_type, bpy.types.Node):
            continue

        # Include geometry nodes and useful shader/function nodes
        if not (
            attr_name.startswith("GeometryNode")
            or attr_name in USEFUL_SHADER_NODES
            or attr_name in USEFUL_FUNCTION_NODES
        ):
            continue

        # Check if it's actually registered
        try:
            if not node_type.is_registered_node_type():
                continue
        except AttributeError:
            continue

        all_nodes.append(node_type)

    # Sort by name for consistent output
    all_nodes.sort(key=lambda x: x.__name__)

    return all_nodes


def generate_file_header() -> str:
    """Generate the header for generated files."""
    return '''"""
Auto-generated Blender Geometry Node classes.

DO NOT EDIT THIS FILE MANUALLY.
This file is generated by molecularnodes/nodes/generator.py

To regenerate: Run generator.py from within Blender

KNOWN LIMITATIONS:
- Dynamic multi-input/output sockets are not yet supported
  (these are the unnamed sockets that appear in the UI for nodes like
  "Evaluate Closure", "Join Geometry", etc. that allow dragging in
  multiple connections)
- TODO: Add support for dynamic socket creation
"""

from __future__ import annotations
from typing import Any
from bpy.types import NodeSocket
from ..builder import NodeBuilder, TreeBuilder

# Type aliases for node inputs
LINKABLE = "NodeSocket | NodeBuilder | Any"
TYPE_INPUT_VECTOR = "tuple[float, float, float] | NodeSocket | NodeBuilder | None"
TYPE_INPUT_ROTATION = "tuple[float, float, float, float] | NodeSocket | NodeBuilder | None"
TYPE_INPUT_BOOLEAN = "bool | NodeSocket | NodeBuilder | None"

'''


def generate_all(output_dir: Path | None = None):
    """Generate all node classes and write to files."""
    if output_dir is None:
        # Default to _generated directory next to this file
        output_dir = Path(__file__).parent / "generated"

    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    print(f"Generating node classes to {output_dir}")

    # Get all geometry nodes
    all_nodes = get_all_geometry_nodes()
    print(f"Found {len(all_nodes)} geometry nodes")

    # Introspect all nodes
    node_infos = []
    for node_type in all_nodes:
        node_info = introspect_node(node_type)
        if node_info:
            node_infos.append(node_info)

    print(f"Successfully introspected {len(node_infos)} nodes")

    # Group by category
    by_category = {}
    for node_info in node_infos:
        category = node_info.category.lower()
        if category not in by_category:
            by_category[category] = []
        by_category[category].append(node_info)

    # Generate files by category
    generated_files = []
    for category, nodes in by_category.items():
        filename = f"{category}.py"
        filepath = output_dir / filename

        print(f"Generating {filename} with {len(nodes)} nodes...")

        with open(filepath, "w") as f:
            f.write(generate_file_header())
            f.write("\n\n")

            for node_info in nodes:
                try:
                    class_code, has_dynamic = generate_node_class(node_info)
                    f.write(class_code)
                    f.write("\n\n")
                except Exception as e:
                    print(f"  Error generating {node_info.name}: {e}")

        generated_files.append(filename)

    # Generate __init__.py that exports everything
    init_file = output_dir / "__init__.py"
    with open(init_file, "w") as f:
        f.write('"""Auto-generated geometry node classes."""\n\n')

        # Import from all category files
        for filename in sorted(generated_files):
            module_name = filename.replace(".py", "")
            f.write(f"from .{module_name} import *\n")

    print("\nGeneration complete!")
    print(f"Generated {len(generated_files)} files:")
    for filename in sorted(generated_files):
        print(f"  - {filename}")
    print(f"\nTotal: {len(node_infos)} node classes")


if __name__ == "__main__":
    generate_all()
