import math
import os
import tempfile

os.environ.setdefault("BLENDER_USER_EXTENSIONS", tempfile.mkdtemp())
import bpy

bpy.ops.wm.read_homefile(use_empty=True)
scene = bpy.context.scene
w = bpy.data.worlds.new("w")
scene.world = w
w.use_nodes = False
w.color = (0.5, 0.5, 0.5)


def mk(name, use_nodes):
    mat = bpy.data.materials.new(name)
    bpy.data.materials.create_gpencil_data(mat)
    mat.grease_pencil.show_fill = True
    mat.grease_pencil.fill_color = (1, 0.5, 0, 1)
    mat.grease_pencil.color = (0, 0, 1, 1)
    mat.use_nodes = use_nodes
    return mat


for x, use_nodes in ((-1.5, True), (1.5, False)):
    gpd = bpy.data.grease_pencils.new("gpd")
    ob = bpy.data.objects.new("GP", gpd)
    scene.collection.objects.link(ob)
    gpd.materials.append(mk(f"m{use_nodes}", use_nodes))
    d = gpd.layers.new("L").frames.new(1).drawing
    d.add_strokes([32])
    s = d.strokes[0]
    s.cyclic = True
    for i, p in enumerate(s.points):
        a = 2 * math.pi * i / 32
        p.position = (x + math.cos(a), math.sin(a), 0)
        p.radius = 0.05
cam = bpy.data.objects.new("cam", bpy.data.cameras.new("cam"))
scene.collection.objects.link(cam)
cam.location = (0, 0, 10)
scene.camera = cam
scene.render.engine = "BLENDER_EEVEE"
scene.render.resolution_x, scene.render.resolution_y = 300, 150
scene.render.filepath = os.path.join(os.path.dirname(__file__), "probe5.png")
bpy.ops.render.render(write_still=True)
m = bpy.data.materials["mFalse"].grease_pencil
for p in m.bl_rna.properties:
    if p.identifier != "rna_type":
        print(p.identifier, getattr(m, p.identifier))
