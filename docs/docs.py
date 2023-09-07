import bpy
from quartodoc import MdRenderer
import griffe
import os
import pathlib
import json
import colorsys


folder = pathlib.Path(__file__).resolve().parent

nodes_json = os.path.join(folder, "nodes.json")
output_file = os.path.join(folder, "nodes.qmd")

# load the long-form descriptions of the nodes from the nodes.json file
with open(nodes_json, 'r') as file:
    node_descriptions = json.load(file)

# load the data file which contains the nodes
mn_data_file = os.path.join(folder, "../MolecularNodes/assets/MN_data_file.blend")
bpy.ops.wm.open_mainfile(filepath = mn_data_file)

def rgba2hex(rgba):
    r = int(rgba[0] * 255)
    g = int(rgba[1] * 255)
    b = int(rgba[2] * 255)
    a = int(rgba[3] * 255)
    return f'#{a:02x}{r:02x}{g:02x}{b:02x}'

def rgb_01_to_hex(r, g, b):
    r, g, b = [int(x * 255) for x in colorsys.rgb_to_hsv(r, g, b)]
    return f'#{r:02x}{g:02x}{b:02x}'

nodes = []

# get the relevant nodes to write documentation about. Currentling limiting it to 
# those that start with "MN". The nodes which start with ".MN" are for interal use only.
for node in bpy.data.node_groups:
    if node.name.startswith("MN"):
        nodes.append(node)

accepted_types = ["Float", "Integer", "Boolean", 'Vector', 'Material', 'Geometry', 'Color']

def get_values(sockets):
    parameter_list = []
    for socket in sockets:
        type = socket.bl_label
        subtype = socket.bl_subtype_label
        
        if type not in accepted_types:
            continue
        
        if type == "Float":
            default = round(socket.default_value, 2)
        elif type == "Geometry":
            default = None
        elif type == "Vector":
            default = [round(value, 2) for value in socket.default_value]
        elif type == "Material":
            default = 'MN_atomic_style'
        elif type == 'Color':
            values = [int(val * 255) for val in list(socket.default_value)[0:3]]
            # color = rgba2hex(list(socket.default_value))
            default = "rgb({}, {}, {})".format(*values)
        else:
            default = socket.default_value
        
        if subtype != 'None':
            type = f"{type} ({subtype})"
        
        parameter_list.append(
            griffe.docstrings.dataclasses.DocstringParameter(
                name = socket.name, 
                annotation = type, 
                value = default, 
                description = socket.description
            )
        )
    return parameter_list

objects = []
cat = ''

text = griffe.docstrings.dataclasses.DocstringSectionText
for node in nodes:
    inputs = griffe.docstrings.dataclasses.DocstringSectionParameters(
        get_values(node.inputs)
    )
    outputs = griffe.docstrings.dataclasses.DocstringSectionParameters(
        get_values(node.outputs)
    )
    
    title = node.name.replace('MN_', '').replace('_', ' ').title()
    if title.split()[0] != cat:
        cat = title.split()[0]
        objects.append(
            [griffe.docstrings.dataclasses.DocstringSectionText(title = None, value = f"## {cat}")]
        )
    desc = node_descriptions.get(node.name)
    
    title = text(title = None, value = f"### {title}")
    inputs_title = text(value = "\n#### Inputs")
    outputs_title = text(value = "\n#### Outputs")
    video = griffe.docstrings.dataclasses.DocstringSectionText(title = None, value = f"![](videos/nodes/{node.name}.mp4)")
    # inputs = griffe.docstrings.dataclasses.DocstringSectionParameters(input_list)
    
    longer_desc = None
    if desc:
        longer_desc = griffe.docstrings.dataclasses.DocstringSectionText(title = None, value = desc)
        objects.append([title, longer_desc, video, inputs_title, inputs, outputs_title, outputs])
    else:
        objects.append([title, video, inputs_title, inputs, outputs_title, outputs])

ren = MdRenderer(header_level = 2)

header = """
---
title: Nodes
toc: true
toc-depth: 3
---
"""

with open(output_file, "w") as file:
    file.write(header)
    for doc in objects:
        section = ''
        for sec in doc:
            text = ren.render(sec)
            
            file.write(text)
            file.write("\n")
            file.write("\n")
        file.write("\n")