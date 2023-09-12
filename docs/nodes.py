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
it is linked and appended to the current working file. By default these node groups
remain linked to the original source file, which mean they can't be manually edited 
(this helps to not accidentally break the internals of MN). 

If you wish to play around and tweak the insides of the nodes yourself, you can do so
by following these steps: 

  1. Right click on the node and select `Show/Hide` -> `Node Options`

  2. Click the `Linked Data Block` icon to make a local copy that is specific to just
this file.

![The linked data block icon](images/nodes/linked_data_block.png){width="400px"}

  3. <kbd>Tab</kbd> in to the node to make your edits. Each node group inside the node
group will also need to be made local if you wish to edit their internals as well.

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