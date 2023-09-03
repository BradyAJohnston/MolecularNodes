import bpy
import numpy as np
from . import coll
from . import nodes
from .obj import create_object
from .obj import add_attribute



bpy.types.Scene.MN_import_star_file_path = bpy.props.StringProperty(
    name = 'star_file_path', 
    description = 'File path for the star file to import.', 
    options = {'TEXTEDIT_UPDATE'}, 
    default = '', 
    subtype = 'FILE_PATH', 
    maxlen = 0
    )
bpy.types.Scene.MN_import_star_file_name = bpy.props.StringProperty(
    name = 'star_file_name', 
    description = 'Name of the created object.', 
    options = {'TEXTEDIT_UPDATE'}, 
    default = 'NewStarInstances', 
    subtype = 'NONE', 
    maxlen = 0
    )




def load_star_file(
    file_path, 
    obj_name = 'NewStarInstances', 
    node_tree = True,
    world_scale =  0.01 
    ):
    import starfile
    from eulerangles import ConversionMeta, convert_eulers
    
    star = starfile.read(file_path, always_dict=True)
    
    star_type = None
    # only RELION 3.1 and cisTEM STAR files are currently supported, fail gracefully
    if 'particles' in star and 'optics' in star:
        star_type = 'relion'
    elif "cisTEMAnglePsi" in star[0]:
        star_type = 'cistem'
    else:
        raise ValueError(
        'File is not a valid RELION>=3.1 or cisTEM STAR file, other formats are not currently supported.'
        )
    
    # Get absolute position and orientations    
    if star_type == 'relion':
        df = star['particles'].merge(star['optics'], on='rlnOpticsGroup')

        # get necessary info from dataframes
        # Standard cryoEM starfile don't have rlnCoordinateZ. If this column is not present 
        # Set it to "0"
        if "rlnCoordinateZ" not in df:
            df['rlnCoordinateZ'] = 0
            
        xyz = df[['rlnCoordinateX', 'rlnCoordinateY', 'rlnCoordinateZ']].to_numpy()
        pixel_size = df['rlnImagePixelSize'].to_numpy().reshape((-1, 1))
        xyz = xyz * pixel_size
        shift_column_names = ['rlnOriginXAngst', 'rlnOriginYAngst', 'rlnOriginZAngst']
        if all([col in df.columns for col in shift_column_names]):
            shifts_ang = df[shift_column_names].to_numpy()
            xyz = xyz - shifts_ang 
        euler_angles = df[['rlnAngleRot', 'rlnAngleTilt', 'rlnAnglePsi']].to_numpy()
        image_id = df['rlnMicrographName'].astype('category').cat.codes.to_numpy()
        
    elif star_type == 'cistem':
        df = star[0]
        df['cisTEMZFromDefocus'] = (df['cisTEMDefocus1'] + df['cisTEMDefocus2']) / 2
        df['cisTEMZFromDefocus'] = df['cisTEMZFromDefocus'] - df['cisTEMZFromDefocus'].median()
        xyz = df[['cisTEMOriginalXPosition', 'cisTEMOriginalYPosition', 'cisTEMZFromDefocus']].to_numpy()
        euler_angles = df[['cisTEMAnglePhi', 'cisTEMAngleTheta', 'cisTEMAnglePsi']].to_numpy()
        image_id = df['cisTEMOriginalImageFilename'].astype('category').cat.codes.to_numpy()

    # coerce starfile Euler angles to Blender convention
    
    target_metadata = ConversionMeta(name='output', 
                                    axes='xyz', 
                                    intrinsic=False,
                                    right_handed_rotation=True,
                                    active=True)
    eulers = np.deg2rad(convert_eulers(euler_angles, 
                               source_meta='relion', 
                               target_meta=target_metadata))

    obj = create_object(obj_name, coll.mn(), xyz * world_scale)

    # create the attribute and add the data for the rotations
    add_attribute(obj, 'MOLRotation', eulers, 'FLOAT_VECTOR', 'POINT')

    # create the attribute and add the data for the image id
    add_attribute(obj, 'MOLIMageId', image_id, 'INT', 'POINT')
    
    # create attribute for every column in the STAR file
    for col in df.columns:
        col_type = df[col].dtype    
        # If col_type is numeric directly add
        if np.issubdtype(col_type, np.number):
            add_attribute(obj, col, df[col].to_numpy().reshape(-1), 'FLOAT', 'POINT')
        
        # If col_type is object, convert to category and add integer values
        elif col_type == object:
            codes = df[col].astype('category').cat.codes.to_numpy().reshape(-1)
            add_attribute(obj, col, codes, 'INT', 'POINT')
            # Add the category names as a property to the blender object
            obj[col + '_categories'] = list(df[col].astype('category').cat.categories)
    
    if node_tree:
        nodes.create_starting_nodes_starfile(obj)
    
    return obj


def panel(layout_function, scene):
    col_main = layout_function.column(heading = "", align = False)
    col_main.label(text = "Import Star File")
    row_import = col_main.row()
    row_import.prop(
        bpy.context.scene, 'MN_import_star_file_name', 
        text = 'Name', 
        emboss = True
    )
    col_main.prop(
        bpy.context.scene, 'MN_import_star_file_path', 
        text = '.star File Path', 
        emboss = True
    )
    row_import.operator('mn.import_star_file', text = 'Load', icon = 'FILE_TICK')



class MN_OT_Import_Star_File(bpy.types.Operator):
    bl_idname = "mn.import_star_file"
    bl_label = "Import Star File"
    bl_description = "Will import the given file, setting up the points to instance an object."
    bl_options = {"REGISTER"}

    @classmethod
    def poll(cls, context):
        return True

    def execute(self, context):
        load_star_file(
            file_path = bpy.context.scene.MN_import_star_file_path, 
            obj_name = bpy.context.scene.MN_import_star_file_name, 
            node_tree = True
        )
        return {"FINISHED"}