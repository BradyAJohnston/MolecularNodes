import bpy
from ..blender import nodes

def build_menu(layout, items):
    for item in items:
        # print(item)
        if item == "break":
            layout.separator()
        elif item['label'] == "custom":
            for button in item['values']:
                item['function'](layout, 
                    label = button['label'], 
                    field = button['field'], 
                    prefix = button['prefix'], 
                    property_id = button['property_id']
                    )
        elif item['name'].startswith("mn."):
            layout.operator(item['name'])
        else:
            menu_item_interface(layout, item['label'], item['name'], item['description'].removesuffix('.'))

def menu_item_interface(layout_function, 
                        label, 
                        node_name, 
                        node_description='Add custom MolecularNodes node group.', 
                        node_link = False
                        ):
    op = layout_function.operator('mn.add_custom_node_group', text = label)
    op.node_label = label
    op.node_name = node_name
    op.node_description = node_description
    op.node_link = node_link

def button_custom_color(layout, label, field, prefix, property_id, starting_value = 0):
    op = layout.operator('mn.color_custom', text = label)
    op.field = field
    op.prefix = prefix
    op.node_property = property_id
    op.node_name = label.lower()
    op.starting_value = starting_value
    op.description = f"Choose individual colors for each {label}"

def button_custom_selection(layout, label, field, prefix, property_id, starting_value = 0):
    op = layout.operator('mn.selection_custom', text = label)
    op.field = field
    op.prefix = prefix
    op.node_property = property_id
    op.node_name = label.lower()
    op.starting_value = starting_value
    op.description = f"Create individual selections for each {label}"

