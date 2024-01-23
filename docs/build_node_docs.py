import bpy
from quartodoc import MdRenderer
import molecularnodes as mn
import griffe
import os
import sys
import pathlib

sys.path.insert(0, os.path.abspath('..'))

folder = pathlib.Path(__file__).resolve().parent
file_output_qmd = os.path.join(folder, "nodes/index.qmd")

# open the data file
bpy.ops.wm.open_mainfile(filepath=mn.blender.nodes.MN_DATA_FILE)


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
                name=socket.name,
                annotation=dtype,
                value=default,
                description=socket.description
            )
        )
    return param_list


cat = ''
text = griffe.docstrings.dataclasses.DocstringSectionText
params = griffe.docstrings.dataclasses.DocstringSectionParameters

categories = {}
for category, node_list in mn.ui.node_info.menu_items.items():
    objects = []
    objects.append(
        [text(title=None, value=f"## {mn.blender.nodes.format_node_name(category)}")])

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
            urls = entry.get('video_url')

            inputs = params(get_values(
                mn.blender.nodes.inputs(bpy.data.node_groups[name])))
            outputs = params(get_values(
                mn.blender.nodes.outputs(bpy.data.node_groups[name])))

            title = mn.blender.nodes.format_node_name(entry.get('label'))
            entry_list.append(text(title=None, value=f"### {title}"))
            if desc:
                entry_list.append(text(title=None, value=desc))
            if urls:
                if not isinstance(urls, list):
                    urls = [urls]
                [
                    entry_list.append(
                        text(title=None, value=f"![]({url}.mp4)")
                    ) for url in urls
                ]

            if len(inputs.as_dict()['value']) > 0:
                entry_list.append(text(value="\n#### Inputs"))
                entry_list.append(inputs)
            if len(outputs.as_dict()['value']) > 0:
                entry_list.append(text(value="\n#### Outputs"))
                entry_list.append(outputs)

            objects.append(entry_list)
    categories[category] = objects

ren = MdRenderer(header_level=2)

header = """---
toc: true
toc-depth: 3
fig-align: center
---
"""

for category, object in categories.items():
    with open(os.path.join(folder, f'nodes/{category}.qmd'), 'w') as file:
        file.write(header)
        for doc in object:
            section = ''
            for sec in doc:
                file.write(ren.render(sec))
                file.write("\n\n")
            file.write("\n")
