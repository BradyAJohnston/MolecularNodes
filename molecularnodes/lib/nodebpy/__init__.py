import bpy

_v = bpy.app.version

if _v >= (5, 2, 0):
    from ._520.src.nodebpy import *
else:
    raise ImportError(f"Unsupported Blender version: {_v}")
