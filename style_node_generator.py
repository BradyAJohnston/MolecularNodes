"""
Style Node Generator

This module provides functionality to extract Style node information from MN_data_file_4.4.blend
and generate Python class files with proper type annotations and documentation.

The generator can operate in two modes:
1. Extract node data to JSON for later processing
2. Generate Python classes directly from the blend file
"""

import json
import bpy
from pathlib import Path
from typing import Dict, List, Any, Optional, Union, Tuple
from dataclasses import dataclass
import molecularnodes as mn


@dataclass
class EnumOption:
    """Represents an enum option with its identifier and name."""
    identifier: str
    name: str
    description: str = ""


@dataclass
class NodeInput:
    """Represents a node input socket with its properties."""
    name: str
    type: str
    default_value: Any
    description: str = ""
    min_value: Optional[float] = None
    max_value: Optional[float] = None
    subtype: Optional[str] = None
    enum_options: Optional[List[EnumOption]] = None


@dataclass
class StyleNodeInfo:
    """Information about a Style node extracted from the blend file."""
    name: str
    description: str
    inputs: List[NodeInput]


def get_python_type_annotation(node_input: NodeInput) -> str:
    """Convert Blender socket type to Python type annotation."""
    type_mapping = {
        'VALUE': 'float',
        'INT': 'int',
        'BOOLEAN': 'bool',
        'VECTOR': 'Tuple[float, float, float]',
        'STRING': 'str',
        'RGBA': 'Tuple[float, float, float, float]',
        'GEOMETRY': 'Any',  # Skip geometry inputs
    }
    
    # Handle enum types
    if node_input.type == 'ENUM' and node_input.enum_options:
        enum_name = generate_enum_class_name(node_input.name)
        return enum_name
    
    return type_mapping.get(node_input.type, 'Any')


def get_default_value(socket) -> Any:
    """Extract the default value from a socket, handling different types."""
    # Get socket type from either old or new interface
    socket_type = getattr(socket, 'type', None)
    if socket_type is None:
        socket_type = getattr(socket, 'socket_type', getattr(socket, 'bl_socket_idname', 'UNKNOWN'))
        # Convert to simple type
        if 'Geometry' in socket_type:
            socket_type = 'GEOMETRY'
        elif 'Float' in socket_type or 'Value' in socket_type:
            socket_type = 'VALUE'
        elif 'Int' in socket_type:
            socket_type = 'INT'
        elif 'Bool' in socket_type:
            socket_type = 'BOOLEAN'
        elif 'Vector' in socket_type:
            socket_type = 'VECTOR'
        elif 'String' in socket_type:
            socket_type = 'STRING'
        elif 'Color' in socket_type or 'RGBA' in socket_type:
            socket_type = 'RGBA'
    
    if socket_type == 'GEOMETRY':
        return None
        
    try:
        # Try different attribute names for default value
        default = getattr(socket, 'default_value', None)
        if default is None:
            default = getattr(socket, 'value', None)
        
        if default is not None:
            # Handle vector types
            if socket_type == 'VECTOR':
                return tuple(default)
            elif socket_type == 'RGBA':
                return tuple(default)
            else:
                return default
            
    except (AttributeError, TypeError):
        pass
    
    # Fallback to type defaults
    type_defaults = {
        'VALUE': 0.0,
        'INT': 0,
        'BOOLEAN': False,
        'VECTOR': (0.0, 0.0, 0.0),
        'STRING': '',
        'RGBA': (1.0, 1.0, 1.0, 1.0),
    }
    return type_defaults.get(socket_type, None)


def extract_style_nodes() -> Dict[str, StyleNodeInfo]:
    """
    Extract information about all Style nodes from the currently loaded blend file.
    
    Returns:
        Dictionary mapping node names to StyleNodeInfo objects
    """
    # Ensure we have the data file loaded
    if not str(mn.assets.MN_DATA_FILE) in bpy.data.filepath:
        bpy.ops.wm.open_mainfile(filepath=str(mn.assets.MN_DATA_FILE))
    
    style_nodes = {}
    
    # Find all nodes that start with "Style "
    for node_group_name in bpy.data.node_groups.keys():
        if node_group_name.startswith("Style "):
            node_group = bpy.data.node_groups[node_group_name]
            
            
            # Extract inputs (skip first Geometry socket)
            inputs = []
            # Use interface.items_tree for newer Blender versions
            if hasattr(node_group, 'interface') and hasattr(node_group.interface, 'items_tree'):
                input_sockets = [item for item in node_group.interface.items_tree if hasattr(item, 'in_out') and item.in_out == 'INPUT']
            else:
                input_sockets = node_group.inputs
            
            for input_socket in input_sockets:
                # Handle different attribute names between old and new Blender interface
                socket_type = getattr(input_socket, 'type', None)
                enum_options = None
                
                if socket_type is None:
                    # For newer interface items, check socket_type or bl_socket_idname
                    socket_type = getattr(input_socket, 'socket_type', getattr(input_socket, 'bl_socket_idname', 'UNKNOWN'))
                    # Convert socket_type to simple type name
                    if 'Geometry' in socket_type:
                        socket_type = 'GEOMETRY'
                    elif 'Float' in socket_type or 'Value' in socket_type:
                        socket_type = 'VALUE'
                    elif 'Int' in socket_type:
                        socket_type = 'INT'
                    elif 'Bool' in socket_type:
                        socket_type = 'BOOLEAN'
                    elif 'Vector' in socket_type:
                        socket_type = 'VECTOR'
                    elif 'String' in socket_type:
                        socket_type = 'STRING'
                    elif 'Color' in socket_type or 'RGBA' in socket_type:
                        socket_type = 'RGBA'
                
                # Check for enum items (newer interface)
                if hasattr(input_socket, 'enum_items') and input_socket.enum_items:
                    socket_type = 'ENUM'
                    enum_options = []
                    for item in input_socket.enum_items:
                        enum_options.append(EnumOption(
                            identifier=item.identifier,
                            name=item.name,
                            description=getattr(item, 'description', '')
                        ))
                
                # Check for enum items (older interface)
                elif hasattr(input_socket, 'enum_items_static') and input_socket.enum_items_static:
                    socket_type = 'ENUM'
                    enum_options = []
                    for item in input_socket.enum_items_static:
                        enum_options.append(EnumOption(
                            identifier=item[0],
                            name=item[1],
                            description=item[2] if len(item) > 2 else ''
                        ))
                
                # Check if this is a NodeSocketMenu with different enum access patterns
                elif 'Menu' in str(socket_type):
                    # Try multiple ways to access menu items
                    enum_items_found = None
                    
                    # Method 1: enum_definition.enum_items
                    if hasattr(input_socket, 'enum_definition') and hasattr(input_socket.enum_definition, 'enum_items'):
                        enum_items_found = input_socket.enum_definition.enum_items
                    
                    # Method 2: direct enum_items
                    elif hasattr(input_socket, 'enum_items'):
                        enum_items_found = input_socket.enum_items
                    
                    # Method 3: default_value_items (for menu sockets)
                    elif hasattr(input_socket, 'default_value_items'):
                        enum_items_found = input_socket.default_value_items
                        
                    # Method 4: Try to get from node tree interface
                    elif hasattr(input_socket, 'interface') and hasattr(input_socket.interface, 'enum_items'):
                        enum_items_found = input_socket.interface.enum_items
                    
                    if enum_items_found:
                        socket_type = 'ENUM'
                        enum_options = []
                        for item in enum_items_found:
                            # Handle different item formats
                            if hasattr(item, 'identifier'):
                                enum_options.append(EnumOption(
                                    identifier=item.identifier,
                                    name=getattr(item, 'name', item.identifier),
                                    description=getattr(item, 'description', '')
                                ))
                            elif isinstance(item, (list, tuple)) and len(item) >= 2:
                                enum_options.append(EnumOption(
                                    identifier=item[0],
                                    name=item[1],
                                    description=item[2] if len(item) > 2 else ''
                                ))
                    else:
                        # Fallback: Use knowledge-based approach for known menu types
                        enum_options = get_known_menu_options(input_socket.name, input_socket.description, input_socket.default_value)
                        if enum_options:
                            socket_type = 'ENUM'
                
                if socket_type != 'GEOMETRY':
                    node_input = NodeInput(
                        name=input_socket.name,
                        type=socket_type,
                        default_value=get_default_value(input_socket),
                        description=getattr(input_socket, 'description', ''),
                        min_value=getattr(input_socket, 'min_value', None),
                        max_value=getattr(input_socket, 'max_value', None),
                        subtype=getattr(input_socket, 'subtype', None),
                        enum_options=enum_options
                    )
                    inputs.append(node_input)
            
            style_info = StyleNodeInfo(
                name=node_group_name,
                description=getattr(node_group, 'description', ''),
                inputs=inputs
            )
            
            style_nodes[node_group_name] = style_info
    
    return style_nodes


def save_style_data_to_json(output_path: Path) -> None:
    """
    Extract style node data and save it to a JSON file.
    
    Args:
        output_path: Path where to save the JSON file
    """
    style_nodes = extract_style_nodes()
    
    # Convert to serializable format
    serializable_data = {}
    for name, info in style_nodes.items():
        serializable_data[name] = {
            'name': info.name,
            'description': info.description,
            'inputs': [
                {
                    'name': inp.name,
                    'type': inp.type,
                    'default_value': inp.default_value,
                    'description': inp.description,
                    'min_value': inp.min_value,
                    'max_value': inp.max_value,
                    'subtype': inp.subtype,
                    'enum_options': [
                        {
                            'identifier': opt.identifier,
                            'name': opt.name,
                            'description': opt.description
                        }
                        for opt in inp.enum_options
                    ] if inp.enum_options else None
                }
                for inp in info.inputs
            ]
        }
    
    with open(output_path, 'w') as f:
        json.dump(serializable_data, f, indent=2, default=str)
    
    print(f"Style node data saved to {output_path}")


def generate_class_name(node_name: str) -> str:
    """Convert node name to Python class name (e.g., 'Style Cartoon' -> 'StyleCartoon')."""
    return ''.join(word.capitalize() for word in node_name.split())


def generate_style_attribute_name(input_name: str) -> str:
    """Convert input name to Python attribute name (e.g., 'Bond Radius' -> 'bond_radius')."""
    return input_name.lower().replace(' ', '_').replace('-', '_')


def generate_enum_class_name(input_name: str) -> str:
    """Generate enum class name (e.g., 'Sphere Geometry' -> 'SphereGeometryEnum')."""
    return ''.join(word.capitalize() for word in input_name.split()) + 'Enum'


def get_known_menu_options(socket_name: str, description: str, default_value: str) -> Optional[List[EnumOption]]:
    """
    Fallback method to provide known enum options for common menu sockets.
    This is used when the Blender interface doesn't expose enum options programmatically.
    """
    # Common menu patterns based on socket names and descriptions
    menu_definitions = {
        "Sphere Geometry": [
            EnumOption("Point", "Point", "Point cloud representation"),
            EnumOption("Instance", "Instance", "Instanced mesh spheres"),
            EnumOption("Mesh", "Mesh", "Realized mesh spheres")
        ],
        "Backbone Shape": [
            EnumOption("Cylinder", "Cylinder", "Cylindrical backbone"),
            EnumOption("Rectangle", "Rectangle", "Rectangular backbone")
        ],
        "Base Shape": [
            EnumOption("Rectangle", "Rectangle", "Rectangular base shape"),
            EnumOption("Cylinder", "Cylinder", "Cylindrical base shape"),
            EnumOption("Square", "Square", "Square base shape")
        ],
        "U Component": [
            EnumOption("Factor", "Factor", "Use factor component"),
            EnumOption("Position", "Position", "Use position component")
        ],
        "Separate By": [
            EnumOption("chain_id", "Chain ID", "Separate by chain identifier"),
            EnumOption("res_id", "Residue ID", "Separate by residue identifier"),
            EnumOption("entity_id", "Entity ID", "Separate by entity identifier")
        ],
        "Color Source": [
            EnumOption("Alpha Carbon", "Alpha Carbon", "Color from alpha carbon atoms"),
            EnumOption("Backbone", "Backbone", "Color from backbone atoms"),
            EnumOption("All", "All", "Color from all atoms")
        ]
    }
    
    # Try exact name match first
    if socket_name in menu_definitions:
        return menu_definitions[socket_name]
    
    # Try pattern matching from description
    if description:
        desc_lower = description.lower()
        if "point cloud" in desc_lower and "instances" in desc_lower and "mesh" in desc_lower:
            return menu_definitions.get("Sphere Geometry")
    
    # If we have a default value, try to infer options
    if default_value:
        if default_value in ["Point", "Instance", "Mesh"]:
            return menu_definitions.get("Sphere Geometry")
        elif default_value in ["Cylinder", "Rectangle"]:
            if "backbone" in socket_name.lower():
                return menu_definitions.get("Backbone Shape")
            elif "base" in socket_name.lower():
                return menu_definitions.get("Base Shape")
        elif default_value == "chain_id":
            return menu_definitions.get("Separate By")
    
    return None


def generate_enum_class(node_input: NodeInput, style_name: str) -> str:
    """Generate an enum class for a socket with enum options."""
    if not node_input.enum_options:
        return ""
    
    enum_name = generate_enum_class_name(node_input.name)
    
    # Generate enum members
    enum_members = []
    for option in node_input.enum_options:
        # Convert identifier to valid Python identifier
        member_name = option.identifier.upper().replace(' ', '_').replace('-', '_')
        enum_members.append(f'    {member_name} = "{option.identifier}"')
    
    enum_class = f'''class {enum_name}(str, Enum):
    """{node_input.description if node_input.description else f'Enum for {node_input.name} in {style_name}'}
    
    Options
    -------
{chr(10).join(f"    {option.identifier} : {option.name}" + (f" - {option.description}" if option.description else "") for option in node_input.enum_options)}
    """
{chr(10).join(enum_members)}
'''
    
    return enum_class


def generate_python_class(style_info: StyleNodeInfo) -> str:
    """
    Generate a Python class string for a Style node.
    
    Args:
        style_info: StyleNodeInfo object containing node data
        
    Returns:
        String containing the Python class definition
    """
    class_name = generate_class_name(style_info.name)
    style_identifier = style_info.name.replace("Style ", "").lower().replace(" ", "_")
    
    # Generate enum classes first
    enum_classes = []
    for inp in style_info.inputs:
        if inp.type == 'ENUM' and inp.enum_options:
            enum_classes.append(generate_enum_class(inp, style_info.name))
    
    # Generate Socket dataclass entries
    socket_entries = []
    init_params = []
    init_assignments = []
    
    for inp in style_info.inputs:
        attr_name = generate_style_attribute_name(inp.name)
        blender_name = inp.name
        
        socket_entries.append(f'        Socket(name="{attr_name}", blendername="{blender_name}"),')
        
        # Generate type annotation and default value
        py_type = get_python_type_annotation(inp)
        default_val = repr(inp.default_value)
        
        # Add parameter with description comment if available
        param_line = f"        {attr_name}: {py_type} = {default_val},"
        if inp.description:
            param_line += f"  # {inp.description}"
        init_params.append(param_line)
        
        init_assignments.append(f"        self.{attr_name} = {attr_name}")
    
    # Generate numpy-style docstring
    class_docstring = style_info.description if style_info.description else f'Style class for {style_info.name}'
    
    # Generate parameters section for numpy-style docstring
    param_docs = []
    for inp in style_info.inputs:
        attr_name = generate_style_attribute_name(inp.name)
        py_type = get_python_type_annotation(inp)
        description = inp.description if inp.description else f'Value for {inp.name}'
        
        # Add enum options to description
        if inp.type == 'ENUM' and inp.enum_options:
            options_text = ", ".join([f"'{opt.identifier}'" for opt in inp.enum_options])
            description += f". Options: {options_text}"
        
        param_docs.append(f"    {attr_name} : {py_type}")
        param_docs.append(f"        {description}")
    
    params_section = '\n'.join(param_docs) if param_docs else "    None"
    
    # Generate the class
    class_code = f'''{"".join(enum_classes)}class {class_name}(StyleBase):
    """{class_docstring}
    
    Parameters
    ----------
{params_section}
    """
    style = "{style_identifier}"
    socketdata: SocketInfo = [
{chr(10).join(socket_entries)}
    ]

    def __init__(
        self,
{chr(10).join(init_params)}
    ):
        """{class_docstring}

        Parameters
        ----------
{params_section}
        """
{chr(10).join(init_assignments)}
'''
    
    return class_code


def generate_style_classes_file(output_path: Path) -> None:
    """
    Generate a complete Python file with all Style classes.
    
    Args:
        output_path: Path where to save the Python file
    """
    # Extract style nodes
    style_nodes = extract_style_nodes()
    
    # Generate file header
    header = '''"""
Auto-generated Style Classes

This file contains Style classes automatically generated from the MN_data_file_4.4.blend file.
Each class represents a Style node and provides a Python interface to configure the node parameters.

Generated classes follow the same pattern as the existing styles in molecularnodes.nodes.styles,
using the Socket dataclass system and StyleBase inheritance.
"""

from dataclasses import dataclass
from enum import Enum
from typing import List, Optional, Tuple, Any
from bpy.types import GeometryNodeGroup

__all__ = [
'''

    # Generate __all__ list
    class_names = [f'    "{generate_class_name(name)}",' for name in style_nodes.keys()]
    all_list = '\n'.join(class_names)
    
    # Add StyleBase and Socket definitions
    base_definitions = '''
]


@dataclass
class Socket:
    """Represents a mapping between a class attribute and a Blender node input socket.

    Attributes:
        name: The name of the attribute in the Style class
        blendername: The corresponding name of the input socket in the Blender node
    """

    name: str
    blendername: Optional[str] = None


SocketInfo = List[Socket]


class StyleBase:
    style: str
    socketdata: SocketInfo = []

    def update_style_node(self, node_style: GeometryNodeGroup):
        """Update the Blender node inputs with values from this style's attributes.

        Args:
            node_style: The Blender GeometryNodeGroup to update
        """
        for input in node_style.inputs:
            if input.type != "GEOMETRY":
                for arg in self.socketdata:
                    if input.name == arg.blendername:
                        input.default_value = getattr(self, arg.name)


'''

    # Generate all classes
    classes = []
    for style_info in style_nodes.values():
        classes.append(generate_python_class(style_info))
    
    # Combine everything
    full_file = header + all_list + base_definitions + '\n\n'.join(classes)
    
    # Write to file
    with open(output_path, 'w') as f:
        f.write(full_file)
    
    print(f"Generated {len(style_nodes)} style classes in {output_path}")


def main():
    """Main function to demonstrate usage."""
    # Set up paths
    docs_dir = Path(__file__).parent
    json_output = docs_dir / "style_nodes_data.json"
    py_output = docs_dir / "generated_style_classes.py"
    
    # Generate JSON data
    save_style_data_to_json(json_output)
    
    # Generate Python classes
    generate_style_classes_file(py_output)
    
    print("Style node generation complete!")


if __name__ == "__main__":
    main()