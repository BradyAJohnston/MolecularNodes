import bpy
from quartodoc import MdRenderer
import griffe
import os
import pathlib

folder = pathlib.Path(__file__).resolve().parent
mn_data_file = os.path.join(folder, "../MolecularNodes/assets/MN_data_file.blend")

bpy.ops.wm.open_mainfile(filepath = mn_data_file)

nodes = []

for node in bpy.data.node_groups:
    if node.name.startswith("MN"):
        nodes.append(node)

objects = []
for node in nodes:
    input_list = []
    for input in node.inputs:
        if input.type in ["VALUE", "INT", "BOOLEAN", 'VECTOR']:
            
            if input.type == "VALUE":
                default = round(input.default_value, 2)
            if input.type == "VECTOR":
                default = list(input.default_value)
            else:
                default = input.default_value
            input_list.append(
                griffe.docstrings.dataclasses.DocstringParameter(
                    name = input.name, 
                    annotation = input.type, 
                    value = default, 
                    description = input.description
                )
            )
    title = node.name.replace('MN_', '').replace('_', ' ').title()
    title = griffe.docstrings.dataclasses.DocstringSectionText(title = None, value = f"### {title}")
    video = griffe.docstrings.dataclasses.DocstringSectionText(title = None, value = f"![](videos/{node.name}.mp4)")
    inputs = griffe.docstrings.dataclasses.DocstringSectionParameters(input_list)
    objects.append([title, video, inputs])

ren = MdRenderer(header_level = 2)

with open("docs/nodes/nodes_auto.qmd", "w") as file:
    file.write(
        """
---
title: Node Documentation
toc: true
toc-depth: 3
---
"""
        )
    # new_section = True
    # sections = ['animation', 'style', 'dna', 'select', 'utils']
    for doc in objects:
        section = ''
        for sec in doc:
            text = ren.render(sec)
            
            file.write(text)
            file.write("\n")
            file.write("\n")
        file.write("\n")