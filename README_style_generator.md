# Style Node Generator

This tool automatically generates Python classes for Style nodes from the MN_data_file_4.4.blend file.

## Overview

The Style Node Generator extracts information about all Style nodes (nodes starting with "Style ") from the blend file and can generate either:

1. **JSON data file** - Raw node information for further processing
2. **Python class file** - Complete class definitions with proper type annotations

## Generated Classes

Each generated class follows the established pattern in `molecularnodes.nodes.styles`:

- Inherits from `StyleBase`
- Uses `Socket` dataclass system for input mapping
- Provides proper type annotations for all parameters
- Includes docstrings with parameter descriptions
- Follows Python naming conventions (e.g., "Bond Radius" → `bond_radius`)

## Usage

### From Command Line (in Blender)

```bash
# Run within Blender to generate both JSON and Python files
blender --python generate_style_docs.py

# Generate only JSON data
blender --python generate_style_docs.py -- --json

# Generate only Python classes  
blender --python generate_style_docs.py -- --python

# List available Style nodes
blender --python generate_style_docs.py -- --list-nodes

# Specify output directory
blender --python generate_style_docs.py -- --output-dir ./output
```

### From Python Script

```python
import bpy
from style_node_generator import (
    extract_style_nodes,
    save_style_data_to_json, 
    generate_style_classes_file
)

# Load the data file (if not already loaded)
bpy.ops.wm.open_mainfile(filepath=str(mn.assets.MN_DATA_FILE))

# Extract node information
style_nodes = extract_style_nodes()
print(f"Found {len(style_nodes)} Style nodes")

# Generate JSON data
save_style_data_to_json("style_data.json")

# Generate Python classes
generate_style_classes_file("style_classes.py")
```

## Example Generated Class

```python
class StyleCartoon(StyleBase):
    """Secondary structure cartoons"""
    style = "cartoon"
    socketdata: SocketInfo = [
        Socket(name="quality", blendername="Quality"),
        Socket(name="dssp", blendername="Peptide DSSP"),
        Socket(name="cylinders", blendername="Peptide Cylinders"),
        # ... more inputs
    ]

    def __init__(
        self,
        quality: int = 2,  # Subdivision quality for generated geometry
        dssp: bool = False,  # Enable DSSP secondary structure calculation
        cylinders: bool = False,  # Use cylinders for helices
        # ... more parameters with type annotations and descriptions
    ):
        """Initialize the StyleCartoon style.
        
        Args:
            quality: Subdivision quality for generated geometry
            dssp: Enable DSSP secondary structure calculation
            cylinders: Use cylinders for helices
            # ... parameter documentation
        """
        self.quality = quality
        self.dssp = dssp
        self.cylinders = cylinders
        # ... assignments
```

## Files

- `style_node_generator.py` - Main generator module with all functionality
- `generate_style_docs.py` - Simple command-line interface
- `README_style_generator.md` - This documentation

## Features

- **Automatic Type Detection** - Converts Blender socket types to proper Python type annotations
- **Default Value Extraction** - Preserves original default values from blend file
- **Description Support** - Includes node and input descriptions in generated classes
- **Naming Conventions** - Converts Blender names to Python-compliant attribute names
- **JSON Export** - Option to save raw data for custom processing
- **Flexible Output** - Generate classes individually or as a complete module

## Integration

The generated classes can be used as drop-in replacements or extensions to the existing style system:

```python
from generated_style_classes import StyleCartoon

# Create a style instance
cartoon_style = StyleCartoon(
    quality=3,
    dssp=True,
    thickness=0.8
)

# Apply to a node group
cartoon_style.update_style_node(node_group)
```