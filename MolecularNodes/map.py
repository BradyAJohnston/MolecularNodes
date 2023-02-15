import MolecularNodes as mn
import bpy
import mrcfile
import numpy as np
import einops
import starfile
from scipy.spatial.transform import Rotation as R

particle_star_file = '/Users/brady/Downloads/starfile_example.star'
star = starfile.read(particle_star_file)

df = star['particles'].merge(star['optics'], on='rlnOpticsGroup')

# get necessary info from dataframes
xyz = df[['rlnCoordinateX', 'rlnCoordinateY', 'rlnCoordinateZ']].to_numpy()
shifts_ang = df[['rlnOriginXAngst', 'rlnOriginYAngst', 'rlnOriginZAngst']].to_numpy()
pixel_size = df['rlnImagePixelSize'].to_numpy()
euler_angles = df[['rlnAngleRot', 'rlnAngleTilt', 'rlnAnglePsi']].to_numpy()

# Get absolute position and orientations
particle_positions = xyz# - (shifts_ang / pixel_size)
rotation_matrices = R.from_euler(
    seq='ZYZ', angles=euler_angles, degrees=True
).inv().as_matrix()

eulers = R.from_matrix(rotation_matrices).as_euler('zyx')

obj = mn.create_object('instances', mn.coll_mn(), xyz / 1e3)

attribute = obj.data.attributes.new('rot', 'FLOAT_VECTOR', 'POINT')
attribute.data.foreach_set('vector', eulers.reshape(len(eulers) * 3))
