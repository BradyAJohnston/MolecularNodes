import bpy
from ..blender import nodes


def build_menu(layout, context, items):
    bob = context.active_object

    for item in items:
        if item == "break":
            layout.separator()
        elif item["label"] == "custom":
            for button in item["values"]:
                row = layout.row()
                item["function"](
                    row,
                    label=button["label"],
                    field=button["field"],
                    dtype=button["dtype"],
                    prefix=button["prefix"],
                    property_id=button["property_id"],
                )
                row.enabled = bool(bob.get(button["property_id"]))
        elif item["name"].startswith("mn."):
            layout.operator(item["name"])
        else:
            label = item["label"]
            name = item["name"]
            description = item["description"].split("\n")[0].removesuffix(".")
            menu_item_interface(layout, label=label, name=name, description=description)


def menu_item_interface(
    layout_function,
    label,
    name,
    description="Add custom MolecularNodes node group.",
    node_link=False,
):
    op = layout_function.operator("mn.add_custom_node_group", text=label)
    op.node_label = label
    op.node_name = name
    op.node_description = description
    op.node_link = node_link


def button_custom_iswitch(
    layout, label, field, dtype, prefix, property_id, starting_value=0
):
    op = layout.operator("mn.iswitch_custom", text=label)
    op.field = field
    op.dtype = dtype
    op.prefix = prefix
    op.node_property = property_id
    op.node_name = label
    op.starting_value = starting_value

    if dtype == "RGBA":
        op.description = f"Choose custom colors for {label}"
    elif dtype == "BOOLEAN":
        op.description = f"Choose custom selections for {label}"
    else:
        raise ValueError(f"Data type currently not supported: {dtype}")
