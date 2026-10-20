"""Render a GN-generated Grease Pencil circle on a GP object, with a World set.

The stroke renders; without a World (probe.py) EEVEE renders nothing.
"""

import os
import tempfile

os.environ.setdefault("BLENDER_USER_EXTENSIONS", tempfile.mkdtemp())
import bpy
from nodebpy import TreeBuilder
from nodebpy import geometry as g

bpy.ops.wm.read_homefile(use_empty=True)
scene = bpy.context.scene
mat = bpy.data.materials.new("gp")
bpy.data.materials.create_gpencil_data(mat)

with TreeBuilder.geometry("probe", is_modifier=True) as tree:
    o = tree.outputs.geometry("Geometry")
    c = g.CurveCircle(radius=1.0) >> g.SetCurveRadius(radius=0.05)
    (
        g.CurvesToGreasePencil(curves=g.GeometryToInstance(c))
        >> g.SetMaterial(material=mat)
        >> g.SetGreasePencilColor(color=(0, 0.2, 1, 1), mode="STROKE")
        >> o
    )

ob = bpy.data.objects.new("GP", bpy.data.grease_pencils.new("gpd"))
scene.collection.objects.link(ob)
ob.data.materials.append(mat)
ob.modifiers.new("GN", "NODES").node_group = bpy.data.node_groups["probe"]
ev = ob.evaluated_get(bpy.context.evaluated_depsgraph_get())
print("geometry", ev.evaluated_geometry())
for layer in ev.data.layers:
    print(layer.name, len(layer.current_frame().drawing.strokes))

cam = bpy.data.objects.new("cam", bpy.data.cameras.new("cam"))
scene.collection.objects.link(cam)
cam.location = (0, 0, 10)
scene.camera = cam
scene.world = bpy.data.worlds.new("w")
scene.render.engine = "BLENDER_EEVEE"
scene.render.resolution_x = scene.render.resolution_y = 200
scene.render.filepath = os.path.join(os.path.dirname(__file__), "probe2.png")
bpy.ops.render.render(write_still=True)
