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
from decimal import DefaultContext
from pathlib import Path
from turtle import Vec2D
from typing import Any
import bpy
from bpy.types import DATA_PT_CURVES_attributes, bpy_prop_array
from mathutils import Euler, Vector
from numpy.ma import isin


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
    # Replace spaces, hyphens, and other non-alphanumeric characters with underscores
    normalized = name.lower()
    normalized = "".join(c if c.isalnum() else "_" for c in normalized)

    # Remove consecutive underscores and leading/trailing underscores
    while "__" in normalized:
        normalized = normalized.replace("__", "_")
    normalized = normalized.strip("_")

    # If the name starts with a digit or is purely numeric, prefix it
    if normalized and (normalized[0].isdigit() or normalized.isdigit()):
        normalized = f"input_{normalized}"

    # If the name is empty or only underscores, provide a fallback
    if not normalized or normalized == "_":
        normalized = "input_socket"

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


def get_socket_type_annotation(socket_info: SocketInfo) -> str:
    """Get the Python type annotation for socket properties."""
    # Map Blender socket types to proper bpy.types annotations
    type_map = {
        "NodeSocketGeometry": "NodeSocket",
        "NodeSocketBool": "bpy.types.NodeSocketBool",
        "NodeSocketVector": "bpy.types.NodeSocketVector",
        "NodeSocketRotation": "bpy.types.NodeSocketRotation",
        "NodeSocketFloat": "bpy.types.NodeSocketFloat",
        "NodeSocketInt": "bpy.types.NodeSocketInt",
        "NodeSocketString": "bpy.types.NodeSocketString",
        "NodeSocketColor": "bpy.types.NodeSocketColor",
        "NodeSocketMaterial": "bpy.types.NodeSocketMaterial",
        "NodeSocketImage": "bpy.types.NodeSocketImage",
        "NodeSocketObject": "bpy.types.NodeSocketObject",
        "NodeSocketCollection": "bpy.types.NodeSocketCollection",
    }

    return type_map.get(socket_info.bl_socket_type, "NodeSocket")


def format_python_value(value: Any) -> str:
    """Format a Python value as a string for code generation."""
    if value is None:
        return "None"
    elif isinstance(value, str):
        return repr(value)
    elif isinstance(value, bool):
        return str(value)
    elif isinstance(value, (int, float)):
        return str(value)
    elif hasattr(value, "__iter__") and not isinstance(value, str):
        # Handle tuples, lists, vectors
        try:
            if len(value) == 3:
                return f"({value[0]}, {value[1]}, {value[2]})"
            elif len(value) == 4:
                return f"({value[0]}, {value[1]}, {value[2]}, {value[3]})"
            else:
                return str(tuple(value))
        except (TypeError, AttributeError):
            return "None"
    else:
        # For other types, try to get a reasonable representation
        try:
            return repr(value)
        except Exception:
            return "None"


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
                # try:
                value = socket.default_value
                if isinstance(value, (Euler, Vector, bpy_prop_array)):
                    value = list(value)
                socket_info.default_value = value
                # except (AttributeError, TypeError):
                #     pass

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
                        default=getattr(node, prop.identifier),
                    )
                )
            elif prop.type in ["BOOLEAN", "INT", "FLOAT", "STRING"]:
                properties.append(
                    PropertyInfo(
                        identifier=prop.identifier,
                        name=prop.name,
                        prop_type=prop.type,
                        default=prop.default,
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

        # Handle special cases for better naming
        method_name = method_name.replace("_", "")
        if method_name == "and":
            method_name = "l_and"
        elif method_name == "or":
            method_name = "l_or"
        elif method_name == "not":
            method_name = "l_not"
        else:
            # Add underscore suffix to avoid Python keyword conflicts for others
            method_name = f"{method_name}"

        # Skip invalid method names
        if not method_name.replace("_", "").replace("l", "").isalnum():
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
                input_params.append(
                    f"{param_name}: {type_hint} = {socket.default_value}"
                )
                # Use the same parameter name as in the constructor
                call_params.append(f"{param_name}={param_name}")

        params_str = ",\n        ".join(input_params)
        call_params_str = ", ".join(call_params)

        # Add operation parameter to call
        operation_param = f'{operation_enum.identifier}="{enum_id}"'
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

        if hasattr(socket, "default_value"):
            default = getattr(socket, "default_value", "None")
            if isinstance(default, str):
                default = f'"{default}"'
        else:
            default = None
        init_params.append(f"{param_name}: {type_hint} = {default}")
        establish_links_params.append((param_name, socket))

    # Add properties as parameters
    for prop in node_info.properties:
        param_name = normalize_name(prop.identifier)
        match prop.prop_type:
            case "ENUM":
                if prop.enum_items:
                    # Create type literal for enum
                    enum_values = [item[0] for item in prop.enum_items]
                    enum_type = (
                        f"Literal[{', '.join(repr(val) for val in enum_values)}]"
                    )
                    init_params.append(f'{param_name}: {enum_type} = "{prop.default}"')
            case "BOOLEAN":
                init_params.append(f"{param_name}: bool  = {prop.default}")
            case "INT":
                init_params.append(f"{param_name}: int  = {prop.default}")
            case "FLOAT":
                init_params.append(f"{param_name}: float  = {prop.default}")
            case "STRING":
                init_params.append(f'{param_name}: str  = "{prop.default}"')
            case _:
                init_params.append(f"{param_name}: Any | None = None")

    # Format init signature
    if len(init_params) > 2:  # If more than just self
        init_signature = (
            "(\n        "
            + ",\n        ".join(init_params)
            + ",\n        **kwargs\n    )"
        )
    else:
        init_signature = "(" + ", ".join(init_params) + ", **kwargs)"

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
            establish_call = f"""        key_args = {{
            {", ".join(link_mappings)}
        }}
        key_args.update(kwargs)
        self._establish_links(**key_args)"""
    else:
        establish_call = "        self._establish_links(**kwargs)"

    # Build property setting calls
    property_calls = []
    for prop in node_info.properties:
        param_name = normalize_name(prop.identifier)
        property_calls.append(f"""        if {param_name} is not None:
            self.node.{prop.identifier} = {param_name}""")

    property_setting = "\n".join(property_calls) if property_calls else ""

    # Generate input properties
    input_properties = []
    used_input_names = set()
    for socket in node_info.inputs:
        if not socket.identifier or socket.identifier.strip() == "":
            continue

        prop_name = f"i_{normalize_name(socket.name)}"
        if prop_name in used_input_names:
            prop_name = f"i_{normalize_name(socket.identifier)}"

        if prop_name in used_input_names:
            continue

        used_input_names.add(prop_name)
        socket_type_annotation = get_socket_type_annotation(socket)

        input_properties.append(f'''
    @property
    def {prop_name}(self) -> {socket_type_annotation}:
        """Input socket: {socket.name}"""
        return self._input("{socket.identifier}")''')

    # Generate output properties
    output_properties = []
    used_output_names = set()
    for socket in node_info.outputs:
        if not socket.identifier or socket.identifier.strip() == "":
            continue

        prop_name = f"o_{normalize_name(socket.name)}"
        if prop_name in used_output_names:
            prop_name = f"o_{normalize_name(socket.identifier)}"

        if prop_name in used_output_names:
            continue

        used_output_names.add(prop_name)
        socket_type_annotation = get_socket_type_annotation(socket)

        output_properties.append(f'''
    @property
    def {prop_name}(self) -> {socket_type_annotation}:
        """Output socket: {socket.name}"""
        return self._output("{socket.identifier}")''')

    # Generate property accessors for node properties
    property_accessors = []
    for prop in node_info.properties:
        prop_name = normalize_name(prop.identifier)
        if prop.prop_type == "ENUM" and prop.enum_items:
            # Create type literal for enum
            enum_values = [item[0] for item in prop.enum_items]
            enum_type = f"Literal[{', '.join(repr(val) for val in enum_values)}]"
            property_accessors.append(f"""
    @property
    def {prop_name}(self) -> {enum_type}:
        return self.node.{prop.identifier}

    @{prop_name}.setter
    def {prop_name}(self, value: {enum_type}):
        self.node.{prop.identifier} = value""")
        elif prop.prop_type == "BOOLEAN":
            property_accessors.append(f"""
    @property
    def {prop_name}(self) -> bool:
        return self.node.{prop.identifier}

    @{prop_name}.setter
    def {prop_name}(self, value: bool):
        self.node.{prop.identifier} = value""")

    # Generate enum convenience methods
    enum_methods = generate_enum_class_methods(node_info)

    # Add node type annotation
    node_type_annotation = (
        f"bpy.types.{node_info.bl_idname}"
        if node_info.bl_idname.startswith(("Geometry", "Function", "Shader"))
        else "bpy.types.Node"
    )

    # Build class
    class_code = f'''
class {class_name}(NodeBuilder):
    """{node_info.description}"""

    name = "{node_info.bl_idname}"
    node: {node_type_annotation}

    def __init__{init_signature}:
        super().__init__()
{establish_call}
{property_setting}
{enum_methods}
{"".join(input_properties)}
{"".join(output_properties)}
{"".join(property_accessors)}
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
import bpy
from typing import Any
from typing_extensions import Literal
from ..builder import NodeBuilder, NodeSocket
from . import types
from .types import LINKABLE, TYPE_INPUT_BOOLEAN, TYPE_INPUT_VECTOR

'''


def get_manually_specified_nodes() -> set[str]:
    """Get the list of manually specified node bl_idnames."""
    manually_specified_nodes = {
        "FunctionNodeRandomValue",
        "ShaderNodeSeparateXYZ",
        "ShaderNodeCombineXYZ",
        "ShaderNodeMix",
        "ShaderNodeMath",
        "FunctionNodeBooleanMath",
    }
    return manually_specified_nodes


def generate_all(output_dir: Path | None = None):
    """Generate all node classes and write to files."""
    if output_dir is None:
        # Default to _generated directory next to this file
        output_dir = Path(__file__).parent / "generated"

    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    print(f"Generating node classes to {output_dir}")

    # Get manually specified nodes to skip
    manually_specified = get_manually_specified_nodes()

    # Get all geometry nodes
    all_nodes = get_all_geometry_nodes()
    print(f"Found {len(all_nodes)} geometry nodes")

    # Introspect all nodes (excluding manually specified ones)
    node_infos = []
    skipped_count = 0
    for node_type in all_nodes:
        if node_type.__name__ in manually_specified:
            print(f"  Skipping manually specified node: {node_type.__name__}")
            skipped_count += 1
            continue

        node_info = introspect_node(node_type)
        if node_info:
            node_infos.append(node_info)

    print(f"Successfully introspected {len(node_infos)} nodes")
    print(f"Skipped {skipped_count} manually specified nodes")

    # Group by category
    by_category = {}
    for node_info in node_infos:
        category = node_info.category.lower()
        if category not in by_category:
            by_category[category] = []
        by_category[category].append(node_info)

    # Generate files by category
    generated_files = []
    NODES_TO_SKIP = ["Closure", "Simulation", "Repeat", "IndexSwitch", "MenuSwitch"]
    for category, nodes in by_category.items():
        filename = f"{category}.py"
        filepath = output_dir / filename

        print(f"Generating {filename} with {len(nodes)} nodes...")

        with open(filepath, "w") as f:
            f.write(generate_file_header())
            f.write("\n\n")

            for node_info in nodes:
                if any([n in node_info.name for n in NODES_TO_SKIP]):
                    continue
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

        # Import manually specified nodes first
        f.write("# Import manually specified nodes\n")
        f.write("from .manually_specified import *\n\n")

        # Import from all category files
        f.write("# Import auto-generated nodes\n")
        for filename in sorted(generated_files):
            module_name = filename.replace(".py", "")
            if module_name != "manually_specified":  # Don't import twice
                f.write(f"from .{module_name} import *\n")

    print("\nGeneration complete!")
    print(f"Generated {len(generated_files)} files:")
    for filename in sorted(generated_files):
        print(f"  - {filename}")
    print(f"\nTotal: {len(node_infos)} node classes")
    print(f"Plus {len(manually_specified)} manually specified nodes")


if __name__ == "__main__":
    generate_all()
