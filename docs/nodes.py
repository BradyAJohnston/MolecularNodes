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
mn_data_file = os.path.join(folder, "../molecularnodes/assets/MN_data_file.blend")
bpy.ops.wm.open_mainfile(filepath = mn_data_file)


# get the relevant nodes to write documentation about. Currentling limiting it to 
# those that start with "MN". The nodes which start with ".MN" are for interal use only.
nodes = []
for node in bpy.data.node_groups:
    if node.name.startswith("MN"):
        nodes.append(node)

accepted_types = ["Float", "Integer", "Boolean", 'Vector', 'Material', 'Geometry', 
                  'Color', 'Collection', 'Object']

def get_values(sockets):
    parameter_list = []
    for socket in sockets:
        type = socket.bl_label
        subtype = socket.bl_subtype_label
        
        if type not in accepted_types:
            continue
        
        if type == "Float":
            default = round(socket.default_value, 2)
        elif type == "Geometry" or type == 'Collection' or type == 'Object':
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
    if node.name.startswith("MN_dna_"):
        continue
    entry_list = []
    name = node.name
    info = node_descriptions.get(name, {'desc': None, 'url': None})
    
    desc = info.get('description')
    url  = info.get('video_url')
    
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
    if url:
        entry_list.append(text(title = None, value = f"![]({url}.mp4)"))

    
    if len(inputs.as_dict()['value']) > 0:
        entry_list.append(text(value = "\n#### Inputs"))
        entry_list.append(inputs)
    if len(outputs.as_dict()['value']) > 0:
        entry_list.append(text(value = "\n#### Outputs"))
        entry_list.append(outputs)
    
    objects.append(entry_list)

ren = MdRenderer(header_level = 2)

header = """
---
title: Nodes
toc: true
toc-depth: 3
fig-align: center
---

The nodes that you work with inside of Molecular Nodes are pre-defined in a 
`MN_data_file.blend` file that is included with the add-on. When a node is added, 
it is appended to the current file to be used by you. If you wish to change how a node
behaves, you can <kbd>Tab</kbd> into the node group to change the internals, and 
<kbd>Ctrl</kbd> + <kbd>Tab</kbd> to leave the node group.

The nodes are shared between node trees, so if you change the internals of the 'Cartoon'
node then this will change how it works for every node tree in your `.blend` file. If
you wish to create changes that are only for a single node tree, you can create a copy 
of it first before making changes by doing the following:

 1. Right click on the node and select `Show / Hide` -> `Node Options`
 2. Click the `Make Copy` icon to the right of the node group's name
 3. It will create a duplicate node group called `MN_style_cartoon.001` for example, 
 which you can now freely edit without changing the other `MN_style_cartoon` node 
 groups.

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