import bpy
import pyopenvdb as vdb
import numpy as np
import mrcfile
import tempfile
import os

def read_map(file):
    """Reads an MRC file and converts it into a pyopenvdb FloatGrid object.

    Args:
        file (str): Path to the MRC file.

    Returns:
        pyopenvdb.FloatGrid: A pyopenvdb FloatGrid object.
    """
    volume = mrcfile.read(file)
    grid = vdb.FloatGrid()
    grid.copyFromArray(volume)
    grid.gridClass = vdb.GridClass.FOG_VOLUME
    grid.name = 'density'
    return grid

def add_density(grid, scale=[1, 1, 1], rotation=[0, 0, 0]):
    """Adds density to the Blender scene using the provided pyopenvdb FloatGrid object.

    Args:
        grid (pyopenvdb.FloatGrid): A pyopenvdb FloatGrid object containing the density data.
        scale (list[float], optional): A list of 3 floats specifying the scaling factor of the density. Defaults to [1, 1, 1].
        rotation (list[float], optional): A list of 3 floats specifying the rotation of the density in radians. Defaults to [0, 0, 0].

    Returns:
        bpy.types.Object: A Blender object containing the density data.
    """
    tmp = tempfile.NamedTemporaryFile(delete=False)
    try:
        vdb.write(tmp.name, grid)
        tmp.close()
        
        
        file_name = os.path.basename(tmp.name).split()[0]
        bpy.ops.object.volume_import(filepath=tmp.name, files=[], scale=scale, rotation=rotation)
        vol = bpy.context.scene.objects[file_name]
        vol.scale = scale
        bpy.ops
    finally:
        os.remove(tmp.name)
    return vol

def load_volume(file: str, name: str = 'Volume', world_scale: float = 0.01) -> bpy.types.Object:
    """Loads an MRC file as a density object into the Blender scene.

    Args:
        file (str): The path to the MRC file.
        name (str, optional): The name of the volume object. Defaults to 'Volume'.
        world_scale (float, optional): The scaling factor to be applied to the voxel size of the density. 
                                       Defaults to 0.01.

    Returns:
        bpy.types.Object: A Blender object containing the density data.
    """
    # Read the MRC file and extract the voxel size
    grid = read_map(file)
    with mrcfile.open(file) as mrc:
        voxel_size = mrc.voxel_size
    
    # Convert the voxel size to a numpy array
    voxel_size = np.array([voxel_size.x, voxel_size.y, voxel_size.z])
    
    # Rotate and scale the grid
    grid.transform.rotate(np.pi / 2, vdb.Axis(1))
    grid.transform.scale(np.array((-1, 1, 1)) * world_scale * voxel_size)
    
    # Add the density to the Blender scene
    vol = add_density(grid)
    
    # Set the name of the volume object
    vol.name = name
    
    return vol