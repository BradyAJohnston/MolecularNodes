import bpy
import numpy as np

from ..blender import (
    nodes, coll, obj
)

__all__ = ['load']

bpy.types.Scene.MN_import_star_file_path = bpy.props.StringProperty(
    name = 'File', 
    description = 'File path for the `.star` file to import.', 
    subtype = 'FILE_PATH', 
    maxlen = 0
    )
bpy.types.Scene.MN_import_star_file_name = bpy.props.StringProperty(
    name = 'Name', 
    description = 'Name of the created object.', 
    default = 'NewStarInstances', 
    maxlen = 0
    )

def load(
    file_path, 
    name = 'NewStarInstances', 
    node_tree = True,
    world_scale =  0.01 
    ):
    import starfile
    from pandas import DataFrame
    from eulerangles import ConversionMeta, convert_eulers
    
    star = starfile.read(file_path)
    star_type = None
    if isinstance(star, dict):
        star = star[star.keys()[0]]
    
    # only RELION 3.1 and cisTEM STAR files are currently supported, fail gracefully
    if 'particles' in star and 'optics' in star:
        star_type = 'relion'
    elif "cisTEMAnglePsi" in star:
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
        if isinstance(star, DataFrame):
            df = star
        else:
            df = star
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

    ensemble = obj.create_object(xyz * world_scale, collection=coll.mn(), name=name)

    # create the attribute and add the data for the rotations
    obj.add_attribute(ensemble, 'MOLRotation', eulers, 'FLOAT_VECTOR', 'POINT')

    # create the attribute and add the data for the image id
    obj.add_attribute(ensemble, 'MOLIMageId', image_id, 'INT', 'POINT')
    
    # create attribute for every column in the STAR file
    for col in df.columns:
        col_type = df[col].dtype    
        # If col_type is numeric directly add
        if np.issubdtype(col_type, np.number):
            obj.add_attribute(ensemble, col, df[col].to_numpy().reshape(-1), 'FLOAT', 'POINT')
        
        # If col_type is object, convert to category and add integer values
        elif col_type == object:
            codes = df[col].astype('category').cat.codes.to_numpy().reshape(-1)
            obj.add_attribute(ensemble, col, codes, 'INT', 'POINT')
            # Add the category names as a property to the blender object
            ensemble[col + '_categories'] = list(df[col].astype('category').cat.categories)
    
    if node_tree:
        nodes.create_starting_nodes_starfile(ensemble)
    
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
            file_path = scene.MN_import_star_file_path, 
            name = scene.MN_import_star_file_name, 
            node_tree = True
        )
        return {"FINISHED"}


def panel(layout, scene):
    layout.label(text = "Load Star File", icon='FILE_TICK')
    layout.separator()
    row_import = layout.row()
    row_import.prop(scene, 'MN_import_star_file_name')
    layout.prop(scene, 'MN_import_star_file_path')
    row_import.operator('mn.import_star_file')