import bpy

from . import parse

bpy.types.Scene.mol_import_cell_pack_path = bpy.props.StringProperty(
    name='File',
    description='File to import (.cif, .bcif)',
    subtype='FILE_PATH',
    maxlen=0
)
bpy.types.Scene.mol_import_cell_pack_name = bpy.props.StringProperty(
    name='Name',
    description='Name of the created object.',
    default='NewCellPackModel',
    maxlen=0
)


def load(
    file_path,
    name='NewCellPackModel',
    node_setup=True,
    world_scale=0.01,
    fraction: float = 1,
):

    ensemble = parse.CellPack(file_path)
    model = ensemble.create_model(
        name=name,
        node_setup=node_setup,
        world_scale=world_scale,
        fraction=fraction
    )

    return model


class MN_OT_Import_Cell_Pack(bpy.types.Operator):
    bl_idname = "mol.import_cell_pack"
    bl_label = "Load"
    bl_description = ""
    bl_options = {"REGISTER"}

    @classmethod
    def poll(cls, context):
        return True

    def execute(self, context):
        s = context.scene
        load(
            file_path=s.mol_import_cell_pack_path,
            name=s.mol_import_cell_pack_name,
            node_setup=True
        )
        return {"FINISHED"}


def panel(layout, scene):
    layout.label(text="Load CellPack Model", icon='FILE_TICK')
    layout.separator()
    row_import = layout.row()
    row_import.prop(scene, 'mol_import_cell_pack_name')
    layout.prop(scene, 'mol_import_cell_pack_path')
    row_import.operator('mol.import_cell_pack')
