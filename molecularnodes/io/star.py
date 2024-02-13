import bpy
from . import parse

bpy.types.Scene.MN_import_star_file_path = bpy.props.StringProperty(
    name='File',
    description='File path for the `.star` file to import.',
    subtype='FILE_PATH',
    maxlen=0
)
bpy.types.Scene.MN_import_star_file_name = bpy.props.StringProperty(
    name='Name',
    description='Name of the created object.',
    default='NewStarInstances',
    maxlen=0
)


def load(
    file_path,
    name='NewStarInstances',
    node_setup=True,
    world_scale=0.01
):

    ensemble = parse.StarFile(file_path)
    ensemble.create_model(name=name, node_setup=node_setup,
                          world_scale=world_scale)

    return ensemble


class MN_OT_Import_Star_File(bpy.types.Operator):
    bl_idname = "mn.import_star_file"
    bl_label = "Load"
    bl_description = "Will import the given file, setting up the points to instance an object."
    bl_options = {"REGISTER"}

    @classmethod
    def poll(cls, context):
        return True

    def execute(self, context):
        scene = context.scene
        load(
            file_path=scene.MN_import_star_file_path,
            name=scene.MN_import_star_file_name,
            node_setup=True
        )
        return {"FINISHED"}


def panel(layout, scene):
    layout.label(text="Load Star File", icon='FILE_TICK')
    layout.separator()
    row_import = layout.row()
    row_import.prop(scene, 'MN_import_star_file_name')
    layout.prop(scene, 'MN_import_star_file_path')
    row_import.operator('mn.import_star_file')
