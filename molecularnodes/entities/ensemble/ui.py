import bpy
from .star import StarFile
from pathlib import Path

from .cellpack import CellPack

bpy.types.Scene.MN_import_star_file_path = bpy.props.StringProperty(
    name="File",
    description="File path for the `.star` file to import.",
    subtype="FILE_PATH",
    maxlen=0,
)


def load_starfile(file_path, node_setup=True, world_scale=0.01):
    ensemble = StarFile.from_starfile(file_path)
    ensemble.create_object(
        name=Path(file_path).name, node_setup=node_setup, world_scale=world_scale
    )

    return ensemble


class MN_OT_Import_Star_File(bpy.types.Operator):
    bl_idname = "mn.import_star_file"
    bl_label = "Load"
    bl_description = (
        "Will import the given file, setting up the points to instance an object."
    )
    bl_options = {"REGISTER"}

    @classmethod
    def poll(cls, context):
        return True

    def execute(self, context):
        scene = context.scene
        load_starfile(
            file_path=scene.MN_import_star_file_path,
            node_setup=True,
        )
        return {"FINISHED"}


def panel_starfile(layout, scene):
    layout.label(text="Load Star File", icon="FILE_TICK")
    layout.separator()
    row_import = layout.row()
    row_import.prop(scene, "MN_import_star_file_path")
    row_import.operator("mn.import_star_file")


bpy.types.Scene.mol_import_cell_pack_path = bpy.props.StringProperty(
    name="File",
    description="File to import (.cif, .bcif)",
    subtype="FILE_PATH",
    maxlen=0,
)
bpy.types.Scene.mol_import_cell_pack_name = bpy.props.StringProperty(
    name="Name",
    description="Name of the created object.",
    default="NewCellPackModel",
    maxlen=0,
)


def load_cellpack(
    file_path,
    name="NewCellPackModel",
    node_setup=True,
    world_scale=0.01,
    fraction: float = 1,
):
    ensemble = CellPack(file_path)
    ensemble.create_object(
        name=name, node_setup=node_setup, world_scale=world_scale, fraction=fraction
    )

    return ensemble


class MN_OT_Import_Cell_Pack(bpy.types.Operator):
    bl_idname = "mol.import_cell_pack"
    bl_label = "Load"
    bl_description = ""
    bl_options = {"REGISTER"}

    def execute(self, context):
        s = context.scene
        load_cellpack(
            file_path=s.mol_import_cell_pack_path,
            name=s.mol_import_cell_pack_name,
            node_setup=True,
        )
        return {"FINISHED"}


def panel_cellpack(layout, scene):
    layout.label(text="Load CellPack Model", icon="FILE_TICK")
    layout.separator()
    row_import = layout.row()
    row_import.prop(scene, "mol_import_cell_pack_name")
    layout.prop(scene, "mol_import_cell_pack_path")
    row_import.operator("mol.import_cell_pack")
