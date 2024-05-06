import bpy
from pathlib import Path
import warnings
from . import parse

bpy.types.Scene.MN_import_local_path = bpy.props.StringProperty(
    name='File',
    description='File path of the structure to open',
    options={'TEXTEDIT_UPDATE'},
    subtype='FILE_PATH',
    maxlen=0
)
bpy.types.Scene.MN_import_local_name = bpy.props.StringProperty(
    name='Name',
    description='Name of the molecule on import',
    options={'TEXTEDIT_UPDATE'},
    default='NewMolecule',
    maxlen=0
)


def load(
    file_path,
    name="Name",
    centre='',
    del_solvent=True,
    style='spheres',
    build_assembly=False
):
    from biotite import InvalidFileError

    suffix = Path(file_path).suffix
    parser = {
        '.pdb': parse.PDB,
        '.pdbx': parse.CIF,
        '.cif': parse.CIF,
        '.bcif': parse.BCIF,
        '.mol': parse.SDF,
        '.sdf': parse.SDF
    }

    if suffix not in parser:
        raise ValueError(
            f"Unable to open local file. Format '{suffix}' not supported.")
    try:
        molecule = parser[suffix](file_path)
    except InvalidFileError:
        molecule = parse.cif.OldCIF(file_path)

    molecule.create_model(
        name=name,
        style=style,
        build_assembly=build_assembly,
        centre=centre,
        del_solvent=del_solvent
    )
    return molecule

# operator that calls the function to import the structure from a local file


class MN_OT_Import_Protein_Local(bpy.types.Operator):
    bl_idname = "mn.import_protein_local"
    bl_label = "Load"
    bl_description = "Open a local structure file"
    bl_options = {"REGISTER", "UNDO"}

    def execute(self, context):
        scene = context.scene
        file_path = scene.MN_import_local_path

        style = scene.MN_import_style
        if not scene.MN_import_node_setup:
            style = None

        centre = scene.MN_centre_type if scene.MN_import_centre else ''

        mol = load(
            file_path=file_path,
            name=scene.MN_import_local_name,
            centre=centre,
            del_solvent=scene.MN_import_del_solvent,
            style=style,
            build_assembly=scene.MN_import_build_assembly
        )

        # return the good news!
        bpy.context.view_layer.objects.active = mol.object
        self.report({'INFO'}, message=f"Imported '{file_path}' as {mol.name}")
        return {"FINISHED"}

    def invoke(self, context, event):
        return self.execute(context)


def panel(layout, scene):

    layout.label(text='Load a Local File', icon='FILE_TICK')
    layout.separator()

    row_name = layout.row(align=False)
    row_name.prop(scene, 'MN_import_local_name')
    row_name.operator('mn.import_protein_local')

    row_import = layout.row()
    row_import.prop(scene, 'MN_import_local_path')
    layout.separator()

    layout.label(text='Options', icon='MODIFIER')
    options = layout.column(align=True)

    row = options.row()
    row.prop(scene, 'MN_import_node_setup', text='')
    col = row.column()
    col.prop(scene, 'MN_import_style')
    col.enabled = scene.MN_import_node_setup

    row_centre = options.row()

    # row_centre.prop(scene, 'MN_import_centre', icon_value=0)
    # row_centre.prop()
    col_centre = row_centre.column()
    col_centre.prop(scene, 'MN_centre_type', text='')
    col_centre.enabled = scene.MN_import_centre
    options.separator()

    grid = options.grid_flow()
    grid.prop(scene, 'MN_import_build_assembly')
    grid.prop(scene, 'MN_import_del_solvent', icon_value=0)
