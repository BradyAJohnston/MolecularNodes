import bpy
import pyopenvdb as vdb
import mrcfile

def read_map(file):
    volume = mrcfile.read(file)
    grid = vdb.FloatGrid()
    grid.copyFromArray(volume)
    grid.gridClass = vdb.GridClass.FOG_VOLUME
    grid.name = 'density'
    return grid

def add_density(grid):
    import os
    import tempfile

    tmp = tempfile.NamedTemporaryFile(delete=False)
    try:
        vdb.write(tmp.name, grid)
        tmp.close()
        vol = bpy.ops.object.volume_import(filepath=tmp.name, files = [])
    finally:
        os.remove(tmp.name)
    return vol

def load_volume(file):
    grid = read_map(file)
    vol = add_density(grid)
    return vol