import bpy
from .star import StarFile
from pathlib import Path

from .cellpack import CellPack
from bpy.props import StringProperty, BoolProperty


def load_starfile(file_path, node_setup=True, world_scale=0.01):
    ensemble = StarFile.from_starfile(file_path)
    ensemble.create_object(
        name=Path(file_path).name, node_setup=node_setup, world_scale=world_scale
    )

    return ensemble


class ImportEnsemble(bpy.types.Operator):
    filepath: StringProperty(  # type: ignore
        name="File",
        description="File path for the `.star` file to import.",
        subtype="FILE_PATH",
        maxlen=0,
    )
    node_setup: BoolProperty(  # type: ignore
        name="Setup Nodes",
        default=True,
        description="Create and set up a Geometry Nodes tree on import",
    )


class MN_OT_Import_Star_File(ImportEnsemble):
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
        load_starfile(
            file_path=self.filepath,
            node_setup=self.node_setup,
        )
        return {"FINISHED"}


def panel_starfile(layout, scene):
    layout.label(text="Load Star File", icon="FILE_TICK")
    layout.separator()
    row_import = layout.row()
    row_import.prop(scene.mn, "import_star_file_path")
    op = row_import.operator("mn.import_star_file")
    op.filepath = scene.mn.import_star_file_path
    op.node_setup = scene.mn.import_node_setup


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


class MN_OT_Import_Cell_Pack(ImportEnsemble):
    bl_idname = "mol.import_cell_pack"
    bl_label = "Load"
    bl_description = ""
    bl_options = {"REGISTER"}

    def execute(self, context):
        load_cellpack(
            file_path=self.filepath,
            name=Path(self.filepath).name,
            node_setup=self.node_setup,
        )
        return {"FINISHED"}


def panel_cellpack(layout, scene):
    layout.label(text="Load CellPack Model", icon="FILE_TICK")
    layout.separator()
    row_import = layout.row()
    row_import.prop(scene, "mol_import_cell_pack_name")
    layout.prop(scene, "mol_import_cell_pack_path")
    op = row_import.operator("mol.import_cell_pack")
    op.filepath = scene.mol_import_cell_pack_path
    op.node_setup = scene.mol_import_node_setup
