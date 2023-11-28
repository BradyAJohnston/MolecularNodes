import bpy
from quartodoc import MdRenderer
import molecularnodes as mn
import griffe
import os
import sys
import pathlib

sys.path.insert(0, os.path.abspath('..'))

folder = pathlib.Path(__file__).resolve().parent
file_output_qmd = os.path.join(folder, "nodes.qmd")

# open the data file
bpy.ops.wm.open_mainfile(filepath = nodes.MN_DATA_FILE)

def col_to_rgb_str(colors):
    values = [int(val * 255) for val in list(colors)]
    return "rgb({}, {}, {})".format(*values[0:3])

def get_values(sockets):
    param_list = []
    
    for socket in sockets.values():
        dtype = socket.bl_socket_idname.replace("NodeSocket", "")
        default = None
        if dtype == "Float":
            default = round(socket.default_value, 2)
        elif dtype in ['Geometry', 'Collection', 'Object']:
            default = None
        elif dtype == "Vector":
            default = [round(x, 2) for x in socket.default_value]
        elif dtype == "Material":
            default = '`MN_atomic_material`'
        elif dtype == "Color":
            default = col_to_rgb_str(socket.default_value)
        else:
            default = socket.default_value
        
        param_list.append(
            griffe.docstrings.dataclasses.DocstringParameter(
                name = socket.name, 
                annotation = dtype, 
                value = default, 
                description = socket.description
            )
        )
    return param_list

objects = []
cat = ''
text = griffe.docstrings.dataclasses.DocstringSectionText
params = griffe.docstrings.dataclasses.DocstringSectionParameters


for category, node_list in mn.ui.node_info.menu_items.items():
    objects.append([text(title=None, value=f"## {category.title()}")])
    
    
    for item in node_list:
        if isinstance(item, str):
            continue
        
        iter_list = [item]
        
        if item['label'] == "custom":
            iter_list = item['values']
        
        
        for entry in iter_list:
            name = entry['name']
            if name.startswith("mn."):
                name = entry['backup']
            
            entry_list = []
            desc = entry.get('description')
            url = entry.get('video_url')
            
            
            inputs  = params(get_values(mn.blender.nodes.inputs(bpy.data.node_groups[name])))
            outputs = params(get_values(mn.blender.nodes.outputs(bpy.data.node_groups[name])))
            
            title = mn.blender.nodes.format_node_name(entry.get('label'))
            entry_list.append(text(title=None, value=f"### {title}"))
            if desc:
                entry_list.append(text(title=None, value=desc))
            if desc:
                entry_list.append(text(title=None, value=f"![]({url}.mp4)"))
            
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