import bpy
import pyopenvdb as vdb
import numpy as np
import mrcfile

def read_map(file):
    volume = mrcfile.read(file)
    grid = vdb.FloatGrid()
    grid.copyFromArray(volume)
    grid.gridClass = vdb.GridClass.FOG_VOLUME
    grid.name = 'density'
    return grid

def add_density(grid, scale = [1,1,1], rotation = [0,0,0]):
    import os
    import tempfile

    tmp = tempfile.NamedTemporaryFile(delete=False)
    try:
        vdb.write(tmp.name, grid)
        tmp.close()
        vol = bpy.ops.object.volume_import(filepath=tmp.name, files = [], scale = [0.001, 0.001, 0.001], rotation=rotation)
    finally:
        os.remove(tmp.name)
    return vol

def load_volume(file, world_scale = 0.001):
    grid = read_map(file)
    with mrcfile.open(file) as mrc:
        voxel_size = mrc.voxel_size 
    
    voxel_size = np.array([voxel_size.x, voxel_size.y, -voxel_size.z])
    
    scale_adjust = np.array([1, 1, -1])
    rotation_adjust = (0, -np.pi / 2, 0)
    scale = voxel_size * world_scale #* scale_adjust
    
    vol = add_density(grid, scale = scale, rotation=rotation_adjust)
    return vol