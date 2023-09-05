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

node_1 = nodes[2]
objects = []
for node in nodes:
    input_list = []
    # input_list.append(
    #     griffe.docstrings.dataclasses.DocstringSectionText(title = None, value = node.name)
    # )
    for input in node.inputs:
        if input.type == "VALUE":
            input_list.append(
                griffe.docstrings.dataclasses.DocstringParameter(
                    name = input.name, 
                    annotation = input.type, 
                    value = round(input.default_value, 2), 
                    description = input.description
                )
            )
    inputs = griffe.docstrings.dataclasses.DocstringSectionParameters(input_list)
    objects.append(inputs)

ren = MdRenderer(header_level = 2)

[print(ren.render(x)) for x in objects]