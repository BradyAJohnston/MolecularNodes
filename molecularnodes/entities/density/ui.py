import bpy
from pathlib import Path
from .mrc import MRC


def load(
    file_path: str,
    name: str = "NewDensity",
    invert: bool = False,
    setup_nodes: bool = True,
    style: str = "density_surface",
    center: bool = False,
    overwrite: bool = False,
):
    density = MRC(
        file_path=file_path, center=center, invert=invert, overwrite=overwrite
    )
    density.create_object(
        name=Path(file_path).name, setup_nodes=setup_nodes, style=style
    )
    return density


class MN_OT_Import_Map(bpy.types.Operator):
    bl_idname = "mn.import_density"
    bl_label = "Load"
    bl_description = "Import a EM density map into Blender"
    bl_options = {"REGISTER"}

    def execute(self, context):
        scene = context.scene
        load(
            file_path=scene.mn.import_density,
            invert=scene.mn.import_density_invert,
            setup_nodes=scene.mn.import_node_setup,
            style=scene.mn.import_density_style,
            center=scene.mn.import_density_center,
        )
        return {"FINISHED"}


def panel(layout, scene):
    layout.label(text="Load EM Map", icon="FILE_TICK")
    layout.separator()

    row = layout.row()
    row.prop(scene.mn, "import_density")
    row.operator("mn.import_density")

    layout.separator()
    col = layout.column()
    col.alignment = "LEFT"
    col.scale_y = 0.5
    label = f"\
    An intermediate file will be created: {scene.mn.import_density}.vdb\
    Please do not delete this file or the volume will not render.\
    Move the original .map file to change this location.\
    "
    for line in label.strip().split("    "):
        col.label(text=line)

    layout.separator()
    layout.label(text="Options", icon="MODIFIER")

    layout.prop(scene.mn, "import_density_invert")
    layout.prop(scene.mn, "import_density_center")
    row = layout.row()
    row.prop(scene.mn, "import_node_setup", text="")
    col = row.column()
    col.prop(scene.mn, "import_density_style")
    col.enabled = scene.mn.import_node_setup
