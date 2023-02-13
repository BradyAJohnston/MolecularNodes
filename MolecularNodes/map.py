import MolecularNodes as mn
import bpy
import mrcfile
import numpy as np

file = "C:\\Users\\BradyJohnston\\Downloads\\emd_28673.map"

with mrcfile.open(file, header_only = True, permissive = True) as mrc:
    voxel_size = (mrc.header.nz, mrc.header.ny, mrc.header.nx)

volume = mrcfile.open(file)

print(voxel_size)

n_points = prod(np.shape(volume.data))

mn.load.create_object('volume', mn.coll_mn(), locations)