import bpy
from quartodoc import MdRenderer
import griffe
import os
import pathlib
import json

with open('docs/nodes.json', 'r') as file:
    node_descriptions = json.load(file)

folder = pathlib.Path(__file__).resolve().parent
mn_data_file = os.path.join(folder, "../MolecularNodes/assets/MN_data_file.blend")

bpy.ops.wm.open_mainfile(filepath = mn_data_file)

nodes = []

for node in bpy.data.node_groups:
    if node.name.startswith("MN"):
        nodes.append(node)

objects = []
cat = ''
for node in nodes:
    input_list = []
    for input in node.inputs:
        if input.type in ["VALUE", "INT", "BOOLEAN", 'VECTOR', 'MATERIAL', 'GEOMETRY']:
            
            if input.type == "VALUE":
                default = round(input.default_value, 2)
            elif input.type == "GEOMETRY":
                default = None
            elif input.type == "VECTOR":
                default = list(input.default_value)
            elif input.type == "MATERIAL":
                default = 'MN_atomic_style'
            else:
                default = input.default_value
            
            type = input.type
            subtype = input.bl_subtype_label
            if subtype != 'None':
                type = f"{type} ({subtype})"
            
            input_list.append(
                griffe.docstrings.dataclasses.DocstringParameter(
                    name = input.name, 
                    annotation = type, 
                    value = default, 
                    description = input.description
                )
            )
    
    title = node.name.replace('MN_', '').replace('_', ' ').title()
    if title.split()[0] != cat:
        cat = title.split()[0]
        objects.append(
            [griffe.docstrings.dataclasses.DocstringSectionText(title = None, value = f"## {cat}")]
        )
    longer_desc = None
    desc = node_descriptions.get(node.name)
    title = griffe.docstrings.dataclasses.DocstringSectionText(title = None, value = f"### {title}")
    video = griffe.docstrings.dataclasses.DocstringSectionText(title = None, value = f"![](videos/nodes/{node.name}.mp4)")
    inputs = griffe.docstrings.dataclasses.DocstringSectionParameters(input_list)
    if desc:
        longer_desc = griffe.docstrings.dataclasses.DocstringSectionText(title = None, value = desc)
        objects.append([title, longer_desc, video, inputs])
    else:
        objects.append([title, video, inputs])

ren = MdRenderer(header_level = 2)

header = """
---
title: Nodes
toc: true
toc-depth: 3
---
"""

with open("docs/nodes.qmd", "w") as file:
    file.write(header)
    for doc in objects:
        section = ''
        for sec in doc:
            text = ren.render(sec)
            
            file.write(text)
            file.write("\n")
            file.write("\n")
        file.write("\n")