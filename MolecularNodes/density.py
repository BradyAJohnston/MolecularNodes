import bpy
import pyopenvdb as vdb
import mrcfile

em_map_file = "C:\\Users\\BradyJohnston\\Downloads\\emd_28673.map"

volume = mrcfile.read(em_map_file)

grid = vdb.FloatGrid()
grid.copyFromArray(volume)
grid.gridClass = vdb.GridClass.FOG_VOLUME
grid.name = 'density'

vol_file = 'C:\\Users\\BradyJohnston\\Desktop\\volume.vdb'
vdb.write(vol_file, grid)
bpy.ops.object.volume_import(filepath=vol_file, files = [])