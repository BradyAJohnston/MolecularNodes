#!/usr/bin/env python3
"""
Simple script to generate Style node documentation.

Usage:
    python generate_style_docs.py [--json] [--python] [--output-dir DIR]

This script will:
1. Load the MN_data_file_4.4.blend file
2. Extract information from all "Style " nodes
3. Generate either JSON data or Python class files (or both)
"""

import argparse
import sys
from pathlib import Path

# Add the current directory to Python path to import our generator
sys.path.insert(0, str(Path(__file__).parent))

from style_node_generator import (
    extract_style_nodes,
    generate_style_classes_file,
    save_style_data_to_json,
)


def main():
    parser = argparse.ArgumentParser(description="Generate Style node documentation")
    parser.add_argument("--json", action="store_true", help="Generate JSON data file")
    parser.add_argument(
        "--python", action="store_true", help="Generate Python class file"
    )
    parser.add_argument("--output-dir", type=str, default=".", help="Output directory")
    parser.add_argument(
        "--list-nodes", action="store_true", help="Just list available Style nodes"
    )

    args = parser.parse_args()

    # Default to both if neither specified
    if not args.json and not args.python and not args.list_nodes:
        args.json = True
        args.python = True

    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)

    try:
        if args.list_nodes:
            print("Extracting Style nodes...")
            style_nodes = extract_style_nodes()
            print(f"\nFound {len(style_nodes)} Style nodes:")
            for name, info in style_nodes.items():
                print(f"  - {name} ({len(info.inputs)} inputs)")
                if info.description:
                    print(f"    Description: {info.description}")
            return

        if args.json:
            json_path = output_dir / "style_nodes_data.json"
            print(f"Generating JSON data to {json_path}...")
            save_style_data_to_json(json_path)

        if args.python:
            py_path = output_dir / "style.py"
            print(f"Generating Python classes to {py_path}...")
            generate_style_classes_file(py_path)

        print("Generation complete!")

    except Exception as e:
        print(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
