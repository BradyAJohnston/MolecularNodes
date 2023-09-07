import bpy
from quartodoc import MdRenderer
import griffe
import os
import pathlib
import json

folder = pathlib.Path(__file__).resolve().parent

file_nodes_json = os.path.join(folder, "nodes.json")
file_output_qmd = os.path.join(folder, "nodes.qmd")

# load the long-form descriptions of the nodes from the nodes.json file
with open(file_nodes_json, 'r') as file:
    node_descriptions = json.load(file)

# load the data file which contains the nodes
mn_data_file = os.path.join(folder, "../MolecularNodes/assets/MN_data_file.blend")
bpy.ops.wm.open_mainfile(filepath = mn_data_file)


# get the relevant nodes to write documentation about. Currentling limiting it to 
# those that start with "MN". The nodes which start with ".MN" are for interal use only.
nodes = []
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
            # doesn't quite print right at the moment, but this will do for now
            values = [int(val * 255) for val in list(socket.default_value)]
            default = "rgb({}, {}, {})".format(*values[0:3])
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
params = griffe.docstrings.dataclasses.DocstringSectionParameters
for node in nodes:
    entry_list = []
    name = node.name
    desc = node_descriptions.get(name)
    
    inputs = params(get_values(node.inputs))
    outputs = params(get_values(node.outputs))
    
    # Format title for nice printing
    title = name.replace('MN_', '').replace('_', ' ').title()
    if title.split()[0] != cat:
        cat = title.split()[0]
        objects.append([text(title = None, value = f"## {cat}")])
    entry_list.append(text(title = None, value = f"### {title}"))
    
    if desc:
        entry_list.append(text(title = None, value = desc))
    entry_list.append(text(title = None, value = f"![](videos/nodes/{name}.mp4)"))
    entry_list.append(text(value = "\n#### Inputs"))
    entry_list.append(inputs)
    entry_list.append(text(value = "\n#### Outputs"))
    entry_list.append(outputs)
    
    objects.append(entry_list)

ren = MdRenderer(header_level = 2)

header = """
---
title: Nodes
toc: true
toc-depth: 3
---
"""

with open(file_output_qmd, "w") as file:
    file.write(header)
    for doc in objects:
        section = ''
        for sec in doc:
            text = ren.render(sec)
            
            file.write(text)
            file.write("\n")
            file.write("\n")
        file.write("\n")