import bpy
# import pyopenvdb as vdb
import numpy as np
import os

def map_to_grid(file: str, invert: bool = False):
    """Reads an MRC file and converts it into a pyopenvdb FloatGrid object.

    This function reads a file in MRC format, and converts it into a pyopenvdb FloatGrid object,
    which can be used to represent volumetric data in Blender.

    Args:
        file (str): The path to the MRC file.
        invert (bool): Whether to invert the data from the grid, defaulting to False. Some file types
        such as EM tomograms have inverted values, where a high value == low density.

    Returns:
        pyopenvdb.FloatGrid: A pyopenvdb FloatGrid object containing the density data.
    """
    import mrcfile
    import pyopenvdb as vdb

    volume = mrcfile.read(file)
    
    dataType = volume.dtype
    
    # enables different grid types
    
    if dataType == "float32" or dataType == "float64":
        grid = vdb.FloatGrid()
    elif dataType == "int8" or dataType == "int16" or dataType == "int32":
        volume = volume.astype('int32')
        grid = vdb.Int32Grid()
    elif dataType == "int64":
        grid = vdb.Int64Grid()

    if invert:
        volume = np.max(volume) - volume
    
    try:
        grid.copyFromArray(volume)
    except ValueError:
        print(f"Grid data type '{volume.dtype}' is an unsupported type.")
    
    grid.gridClass = vdb.GridClass.FOG_VOLUME
    grid.name = 'density'
    return grid

def path_to_vdb(file: str):
    # Set up file paths
    folder_path = os.path.dirname(file)
    name = os.path.basename(file).split(".")[0]
    file_name = name + '.vdb'
    file_path = os.path.join(folder_path, file_name)
    return file_path
    

def map_to_vdb(file: str, invert: bool = False, world_scale=0.01, overwrite=False) -> str:
    """
    Converts an MRC file to a .vdb file using pyopenvdb.

    Args:
        file (str): The path to the input MRC file.
        invert (bool): Whether to invert the data from the grid, defaulting to False. Some file types
        such as EM tomograms have inverted values, where a high value == low density.
        world_scale (float, optional): The scaling factor to apply to the voxel size of the input file. Defaults to 0.01.
        overwrite (bool, optional): If True, the .vdb file will be overwritten if it already exists. Defaults to False.

    Returns:
        str: The path to the converted .vdb file.
    """
    import mrcfile
    import pyopenvdb as vdb

    file_path = path_to_vdb(file)
    
    # If the map has already been converted to a .vdb and overwrite is False, return that instead
    if os.path.exists(file_path) and not overwrite:
        return file_path

    # Read in the MRC file and convert it to a pyopenvdb grid
    grid = map_to_grid(file, invert = invert)
    
    # Read the voxel size from the MRC file and convert it to a numpy array
    with mrcfile.open(file) as mrc:
        voxel_size = np.array([mrc.voxel_size.x, mrc.voxel_size.y, mrc.voxel_size.z])
    
    # Rotate and scale the grid for import into Blender
    grid.transform.rotate(np.pi / 2, vdb.Axis(1))
    grid.transform.scale(np.array((-1, 1, 1)) * world_scale * voxel_size)
    
    # Write the grid to a .vdb file
    vdb.write(file_path, grid)
    
    # Return the path to the output file
    return file_path

def vdb_to_volume(file: str) -> bpy.types.Object:
    """Imports a VDB file as a Blender volume object.

    Args:
        file (str): Path to the VDB file.

    Returns:
        bpy.types.Object: A Blender object containing the imported volume data.
    """
    # extract name of file for object name
    name = os.path.basename(file).split('.')[0]
    
    # import the volume object
    bpy.ops.object.volume_import(
        filepath = file, 
        files = [], 
        scale = [1, 1, 1], 
        rotation = [0, 0, 0]
    )
    
    # get reference to imported object and return
    vol = bpy.context.scene.objects[name]
    return vol


def load(file: str, name: str = None, invert: bool = False, world_scale: float = 0.01) -> bpy.types.Object:
    """
    Loads an MRC file into Blender as a volumetric object.

    Args:
        file (str): Path to the MRC file.
        name (str, optional): If not None, renames the object with the new name.
        invert (bool): Whether to invert the data from the grid, defaulting to False. Some file types
        such as EM tomograms have inverted values, where a high value == low density.
        world_scale (float, optional): Scale of the object in the world. Defaults to 0.01.

    Returns:
        bpy.types.Object: The loaded volumetric object.
    """
    # Convert MRC file to VDB format
    vdb_file = map_to_vdb(file, invert = invert, world_scale = world_scale)
    
    # Import VDB file into Blender
    vol_object = vdb_to_volume(vdb_file)
    
    if name:
        # Rename object to specified name
        vol_object.name = name
    
    return vol_object
